#include "PCH.h"
#include <numeric>
#include "rs2_utils.h"
#include "rs2_colorizer.h"
#include "rs2_depth_metrics.h"

#define MAX_BUFFER_U16 0xFFFF

template<typename T>
inline UTexture2D* Get(TUniquePtr<T>& Dtex) { return Dtex.Get() ? Dtex.Get()->GetTextureObject() : nullptr; }

inline void TickTex(FDynamicTexture* Tex) { if (Tex) Tex->Tick_GameThread(); }

ARealSenseInspector::ARealSenseInspector(const FObjectInitializer& ObjectInitializer) : Super(ObjectInitializer)
{
	REALSENSE_TRACE(TEXT("+ARealSenseInspector %p"), this);

	PrimaryActorTick.bCanEverTick = true;

	RootComponent = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));

	PclMesh = CreateDefaultSubobject<URuntimeMeshComponent>(TEXT("RuntimeMeshComponent"));
	PclMesh->SetVisibility(false);
	PclMesh->SetCollisionProfileName(UCollisionProfile::NoCollision_ProfileName);
	PclMesh->Mobility = EComponentMobility::Movable;
	PclMesh->bGenerateOverlapEvents = false;
	PclMesh->SetupAttachment(RootComponent);
}

ARealSenseInspector::~ARealSenseInspector()
{
	REALSENSE_TRACE(TEXT("~ARealSenseInspector %p"), this);

	Stop();
}

bool ARealSenseInspector::Start()
{
	SCOPED_PROFILER;

	try
	{
		FScopeLock Lock(&StateMx);

		REALSENSE_TRACE(TEXT("ARealSenseInspector::Start %p"), this);

		if (RsPipeline.Get() || StartedFlag)
		{
			throw std::runtime_error("Already started");
		}

		PclMesh->SetVisibility(false);

		Context = IRealSensePlugin::Get().GetContext();
		if (!Context)
		{
			throw std::runtime_error("GetContext failed");
		}

		rs2::context_ref RsContext(Context->GetHandle());
		rs2::config RsConfig;

		const bool IsPlaybackMode = (PipelineMode == ERealSensePipelineMode::PlaybackFile);
		if (IsPlaybackMode)
		{
			REALSENSE_TRACE(TEXT("enable_device_from_file: %s"), *CaptureFile);
			RsConfig.enable_device_from_file(TCHAR_TO_ANSI(*CaptureFile));
		}
		else
		{
			if (Context->Devices.Num() == 0)
			{
				Context->QueryDevices();
				if (Context->Devices.Num() == 0)
				{
					throw std::runtime_error("No devices available");
				}
			}

			ActiveDevice = DeviceSerial.IsEmpty() ? Context->GetDeviceById(0) : Context->FindDeviceBySerial(DeviceSerial);
			if (!ActiveDevice)
			{
				throw std::runtime_error("Device not found");
			}

			REALSENSE_TRACE(TEXT("enable_device: SN=%s"), *ActiveDevice->Serial);
			RsConfig.enable_device(std::string(TCHAR_TO_ANSI(*ActiveDevice->Serial)));
		}

		const bool EnableDepthStream = bEnableRawDepth || bEnableColorizedDepth || bEnablePcl;
		if (EnableDepthStream)
		{
			if (!IsPlaybackMode)
			{
				REALSENSE_TRACE(TEXT("enable_stream DEPTH %dx%d @%d"), DepthConfig.Width, DepthConfig.Height, DepthConfig.Rate);

				EnsureProfileSupported(ActiveDevice, ERealSenseStreamType::STREAM_DEPTH, ERealSenseFormatType::FORMAT_Z16, DepthConfig);
				RsConfig.enable_stream(RS2_STREAM_DEPTH, DepthConfig.Width, DepthConfig.Height, RS2_FORMAT_Z16, DepthConfig.Rate);

				if (bAlignDepthToColor)
				{
					RsAlign.Reset(new rs2::align(RS2_STREAM_COLOR));
				}
			}

			if (bEnableRawDepth)
			{
				// PF_G16 = DXGI_FORMAT_R16_UNORM ; A single-component, 16-bit unsigned-normalized-integer format that supports 16 bits for the red channel.
				// PF_R16F = DXGI_FORMAT_R16_FLOAT ; A single-component, 16-bit floating-point format that supports 16 bits for the red channel.
				// PF_R16_UINT = DXGI_FORMAT_R16_UINT ; A single-component, 16-bit unsigned-integer format that supports 16 bits for the red channel.
				DepthRawDtex.Reset(new FDynamicTexture(TEXT("DepthRaw"), DepthConfig.Width, DepthConfig.Height, PF_G16, TextureCompressionSettings::TC_Displacementmap));
			}

			if (bEnableColorizedDepth)
			{
				DepthColorizedDtex.Reset(new FDynamicTexture(TEXT("DepthColorized"), DepthConfig.Width, DepthConfig.Height, PF_R8G8B8A8, TextureCompressionSettings::TC_VectorDisplacementmap));
			}
		}

		if (bEnableColor || bAlignDepthToColor)
		{
			if (!IsPlaybackMode)
			{
				REALSENSE_TRACE(TEXT("enable_stream COLOR %dx%d @%d"), ColorConfig.Width, ColorConfig.Height, ColorConfig.Rate);

				EnsureProfileSupported(ActiveDevice, ERealSenseStreamType::STREAM_COLOR, ERealSenseFormatType::FORMAT_RGBA8, ColorConfig);
				RsConfig.enable_stream(RS2_STREAM_COLOR, ColorConfig.Width, ColorConfig.Height, RS2_FORMAT_RGBA8, ColorConfig.Rate);
			}

			if (bEnableColor)
			{
				ColorDtex.Reset(new FDynamicTexture(TEXT("Color"), ColorConfig.Width, ColorConfig.Height, PF_R8G8B8A8, TextureCompressionSettings::TC_VectorDisplacementmap));
			}
		}

		if (bEnableInfrared)
		{
			if (EnableDepthStream)
			{
				InfraredConfig = DepthConfig;
			}

			if (!IsPlaybackMode)
			{
				REALSENSE_TRACE(TEXT("enable_stream INFRARED %dx%d @%d"), InfraredConfig.Width, InfraredConfig.Height, InfraredConfig.Rate);

				EnsureProfileSupported(ActiveDevice, ERealSenseStreamType::STREAM_INFRARED, ERealSenseFormatType::FORMAT_Y8, InfraredConfig);
				RsConfig.enable_stream(RS2_STREAM_INFRARED, InfraredConfig.Width, InfraredConfig.Height, RS2_FORMAT_Y8, InfraredConfig.Rate);
			}

			// PF_G8 = DXGI_FORMAT_R8_UNORM ; A single-component, 8-bit unsigned-normalized-integer format that supports 8 bits for the red channel.
			InfraredDtex.Reset(new FDynamicTexture(TEXT("Infrared"), InfraredConfig.Width, InfraredConfig.Height, PF_G8, TextureCompressionSettings::TC_Displacementmap));
		}

		if (bEnablePcl)
		{
			RsPointCloud.Reset(new rs2::pointcloud());
			RsPoints.Reset(new rs2::points());
		}
		PclCalculateFlag = bEnablePcl;

		if (PipelineMode == ERealSensePipelineMode::RecordFile)
		{
			REALSENSE_TRACE(TEXT("enable_record_to_file: %s"), *CaptureFile);
			RsConfig.enable_record_to_file(TCHAR_TO_ANSI(*CaptureFile));
		}

		FlushRenderingCommands(); // wait for textures
		DepthRawTexture = Get(DepthRawDtex);
		DepthColorizedTexture = Get(DepthColorizedDtex);
		ColorTexture = Get(ColorDtex);
		InfraredTexture = Get(InfraredDtex);

		REALSENSE_TRACE(TEXT("Start pipeline"));
		RsPipeline.Reset(new rs2::pipeline());
		rs2::pipeline_profile RsProfile = RsPipeline->start(RsConfig);
		auto Device = RsProfile.get_device();
		GetDeviceInfo(&Device);
		StartedFlag = true;

		REALSENSE_TRACE(TEXT("Start worker"));
		Worker.Reset(new FRealSenseInspectorWorker(this));
		FString ThreadName(FString::Printf(TEXT("FRealSenseInspectorWorker_%s"), *FGuid::NewGuid().ToString()));
		Thread.Reset(FRunnableThread::Create(Worker.Get(), *ThreadName, 0, TPri_Normal));
		if (!Thread.Get())
		{
			throw std::runtime_error("CreateThread failed");
		}
	}
	catch (const rs2::error& ex)
	{
		auto what(uestr(ex.what()));
		auto func(uestr(ex.get_failed_function()));
		auto args(uestr(ex.get_failed_args()));

		REALSENSE_ERR(TEXT("ARealSenseInspector::Start exception: %s (FUNC %s ; ARGS %s ; TYPE %d"), *what, *func, *args, (int)ex.get_type());
		StartedFlag = false;
	}
	catch (const std::exception& ex)
	{
		REALSENSE_ERR(TEXT("ARealSenseInspector::Start exception: %s"), ANSI_TO_TCHAR(ex.what()));
		StartedFlag = false;
	}

	if (!StartedFlag)
	{
		Stop();
	}

	return StartedFlag ? true : false;
}

void ARealSenseInspector::Stop()
{
	SCOPED_PROFILER;

	try
	{
		FScopeLock Lock(&StateMx);

		StartedFlag = false;

		if (RsPipeline.Get())
		{
			REALSENSE_TRACE(TEXT("Stop pipeline"));
			try {
				NAMED_PROFILER("rs2::pipeline::stop");
				RsPipeline->stop();
			}
			catch (...) {}
		}

		if (Thread.Get())
		{
			REALSENSE_TRACE(TEXT("JoinRealSenseThread"));
			{
				NAMED_PROFILER("JoinRealSenseThread");
				Thread->WaitForCompletion();
			}

			REALSENSE_TRACE(TEXT("Reset thread"));
			Thread.Reset();
		}

		Worker.Reset();
		RsPipeline.Reset();
		RsAlign.Reset();
		RsPoints.Reset();
		RsPointCloud.Reset();

		REALSENSE_TRACE(TEXT("FlushRenderingCommands"));
		{
			NAMED_PROFILER("FlushRenderingCommands");
			FlushRenderingCommands();
		}

		REALSENSE_TRACE(TEXT("DestroyDynamicTextures"));
		{
			NAMED_PROFILER("DestroyDynamicTextures");
			DepthRawDtex.Reset();
			DepthColorizedDtex.Reset();
			ColorDtex.Reset();
			InfraredDtex.Reset();
		}

		ActiveDevice = nullptr;
		DepthRawTexture = nullptr;
		DepthColorizedTexture = nullptr;
		ColorTexture = nullptr;
		InfraredTexture = nullptr;
	}
	catch (const std::exception& ex)
	{
		REALSENSE_ERR(TEXT("ARealSenseInspector::Stop exception: %s"), ANSI_TO_TCHAR(ex.what()));
	}

	#if defined(PROF_ENABLED)
	Profiler::GetInstance()->LogSummary();
	#endif
}

void ARealSenseInspector::ThreadProc()
{
	REALSENSE_TRACE(TEXT("Enter ARealSenseInspector::ThreadProc %p"), this);

	try
	{
		while (StartedFlag)
		{
			if (bEnablePolling)
			{
				PollFrames();
			}
			else
			{
				WaitFrames();
			}
		}
	}
	catch (const std::exception& ex)
	{
		REALSENSE_ERR(TEXT("ARealSenseInspector::ThreadProc exception: %s"), ANSI_TO_TCHAR(ex.what()));
		StartedFlag = false;
	}

	REALSENSE_TRACE(TEXT("Leave ARealSenseInspector::ThreadProc %p"), this);
}

void ARealSenseInspector::PollFrames()
{
	const double t0 = FPlatformTime::Seconds();

	rs2::frameset Frameset;
	bool GotFrames;

	try
	{
		NAMED_PROFILER("rs2::pipeline::poll_for_frames");
		GotFrames = RsPipeline->poll_for_frames(&Frameset);
	}
	catch (const rs2::error& ex)
	{
		REALSENSE_TRACE(TEXT("poll_for_frames failed: %s"), ANSI_TO_TCHAR(ex.what()));
		GotFrames = false;
	}

	if (GotFrames)
	{
		ProcessFrameset(&Frameset);
	}

	const double t1 = FPlatformTime::Seconds();
	const double dt = t1 - t0;
	if (dt < PollFrameRate)
	{
		FPlatformProcess::Sleep(PollFrameRate - dt);
	}
}

void ARealSenseInspector::WaitFrames()
{
	rs2::frameset Frameset;
	bool GotFrames;

	try
	{
		NAMED_PROFILER("rs2::pipeline::wait_for_frames");
		Frameset = RsPipeline->wait_for_frames((unsigned int)(WaitFrameTimeout * 1000.0f));
		GotFrames = true;
	}
	catch (const rs2::error& ex)
	{
		REALSENSE_TRACE(TEXT("wait_for_frames failed: %s"), ANSI_TO_TCHAR(ex.what()));
		GotFrames = false;
	}

	if (GotFrames)
	{
		ProcessFrameset(&Frameset);
	}
}

void ARealSenseInspector::ProcessFrameset(rs2::frameset* Frameset)
{
	SCOPED_PROFILER;

	FScopedTryLock PclScopedMx;

	if (!(bEnableRawDepth && bAnalyzeDepthQuality))
	{
		ResetDepthQuality();
	}

	if (bEnableRawDepth || bEnableColorizedDepth || bEnablePcl)
	{
		rs2::depth_frame DepthFrame = (bAlignDepthToColor && RsAlign.Get()) ? RsAlign->process(*Frameset).get_depth_frame() : Frameset->get_depth_frame();

		if (bEnableRawDepth && DepthRawDtex.Get())
		{
			NAMED_PROFILER("UpdateDepthRaw");
			DepthRawDtex->EnqueUpdateFrame(DepthFrame);
		}

		if (bEnableRawDepth && bAnalyzeDepthQuality)
		{
			AnalyzeDepthQuality(&DepthFrame);
		}

		if (bEnableColorizedDepth && DepthColorizedDtex.Get())
		{
			NAMED_PROFILER("UpdateDepthColorized");
			auto* Tud = DepthColorizedDtex->AllocBuffer();
			rs2::colorize_depth((rs2::pixel_rgba8*)Tud->Data, DepthFrame, (int)DepthColormap, DepthMin, DepthMax, DepthScale, bEqualizeHistogram);
			DepthColorizedDtex->EnqueUpdateData(Tud);
		}

		if (bEnablePcl && RsPointCloud.Get() && PclCalculateFlag)
		{
			if (PclScopedMx.TryLock(&PointCloudMx))
			{
				NAMED_PROFILER("rs2::pointcloud::calculate");
				*RsPoints = RsPointCloud->calculate(DepthFrame);
			}
		}
	}

	if (bEnableColor)
	{
		auto ColorFrame = Frameset->get_color_frame();

		if (ColorDtex.Get())
		{
			NAMED_PROFILER("UpdateColor");
			ColorDtex->EnqueUpdateFrame(ColorFrame);
		}

		if (PclScopedMx.IsLocked())
		{
			NAMED_PROFILER("rs2::pointcloud::map_to");
			RsPointCloud->map_to(ColorFrame);
		}
	}

	if (PclScopedMx.IsLocked())
	{
		PclScopedMx.Unlock();
		PclFramesetId = FramesetId;
		PclCalculateFlag = false;
		PclReadyFlag = true;
	}

	if (bEnableInfrared && InfraredDtex.Get())
	{
		NAMED_PROFILER("UpdateInfrared");
		InfraredDtex->EnqueUpdateFrame(Frameset->get_infrared_frame());
	}

	//REALSENSE_TRACE(TEXT("Frameset %d"), FramesetId);
	FramesetId++;
}

void ARealSenseInspector::Tick(float DeltaSeconds)
{
	SCOPED_PROFILER;

	Super::Tick(DeltaSeconds);

	{
		NAMED_PROFILER("TickTextures");
		TickTex(DepthRawDtex.Get());
		TickTex(DepthColorizedDtex.Get());
		TickTex(ColorDtex.Get());
		TickTex(InfraredDtex.Get());
	}

	if (bEnablePcl && PclMesh && StartedFlag)
	{
		PclRenderAccum += DeltaSeconds;
		if (PclRenderAccum >= PclRenderRate)
		{
			PclRenderAccum = 0;
			PclCalculateFlag = true;
		}

		if (PclReadyFlag)
		{
			FScopedTryLock Mx;
			if (Mx.TryLock(&PointCloudMx))
			{
				UpdatePointCloud();
				PclReadyFlag = false;
			}
		}
	}
	else if (PclMesh)
	{
		PclMesh->SetVisibility(false);
	}
}

void ARealSenseInspector::AnalyzeDepthQuality(rs2::depth_frame* Frame)
{
	SCOPED_PROFILER;

	auto profile = Frame->get_profile();
	auto stream_type = profile.stream_type();

	auto depth_profile = profile.as<rs2::video_stream_profile>();
	rs2_intrinsics intrin = depth_profile.get_intrinsics();
	rs2::region_of_interest roi = { 
		int(intrin.width * (0.5f - 0.5f * AnalyzeROI)),
		int(intrin.height * (0.5f - 0.5f * AnalyzeROI)),
		int(intrin.width * (0.5f + 0.5f * AnalyzeROI)),
		int(intrin.height * (0.5f + 0.5f * AnalyzeROI))
	};

	std::vector<rs2::depth_quality::single_metric_data> sample;

	auto metrics = rs2::depth_quality::analyze_depth_image(
		*Frame, DepthScale, StereoBaseline, &intrin, roi, AnalyzeGroundTruth, DepthQuality.bPlaneFit, sample, false,
		[&](
			const std::vector<rs2::float3>& points,
			const rs2::plane p,
			const rs2::region_of_interest roi,
			const float baseline_mm,
			const float focal_length_pixels,
			const int ground_truth_mm,
			const bool plane_fit,
			const float plane_fit_to_ground_truth_mm,
			const float distance_mm,
			bool record,
			std::vector<rs2::depth_quality::single_metric_data>& samples)
		{
			float TO_METERS = DepthScale;
			static const float TO_MM = 1000.f;
			static const float TO_PERCENT = 100.f;

			// Calculate fill rate relative to the ROI
			auto fill_rate = points.size() / float((roi.max_x - roi.min_x)*(roi.max_y - roi.min_y)) * TO_PERCENT;
			//fill->add_value(fill_rate);
			//if (record) samples.push_back({ fill->get_name(),  fill_rate });
			DepthQuality.FillRate = fill_rate;

			if (!plane_fit)
			{
				DepthQuality.PlaneFitError = -1;
				DepthQuality.SubpixelError = -1;
				DepthQuality.ZAccuracy = -1;
				return;
			}

			const float bf_factor = baseline_mm * focal_length_pixels * TO_METERS; // also convert point units from mm to meter

			std::vector<rs2::float3> points_set = points;
			std::vector<float> distances;
			std::vector<float> disparities;
			std::vector<float> gt_errors;

			// Reserve memory for the data
			distances.reserve(points.size());
			disparities.reserve(points.size());
			if (ground_truth_mm) gt_errors.reserve(points.size());

			// Remove outliers [below 0.5% and above 99.5%)
			std::sort(points_set.begin(), points_set.end(), [](const rs2::float3& a, const rs2::float3& b) { return a.z < b.z; });
			size_t outliers = points_set.size() / 200;
			points_set.erase(points_set.begin(), points_set.begin() + outliers); // crop min 0.5% of the dataset
			points_set.resize(points_set.size() - outliers); // crop max 0.5% of the dataset

															 // Convert Z values into Depth values by aligning the Fitted plane with the Ground Truth (GT) plane
															 // Calculate distance and disparity of Z values to the fitted plane.
															 // Use the rotated plane fit to calculate GT errors
			for (auto point : points_set)
			{
				// Find distance from point to the reconstructed plane
				auto dist2plane = p.a*point.x + p.b*point.y + p.c*point.z + p.d;
				// Project the point to plane in 3D and find distance to the intersection point
				rs2::float3 plane_intersect = { 
					float(point.x - dist2plane*p.a),
					float(point.y - dist2plane*p.b),
					float(point.z - dist2plane*p.c) };

				// Store distance, disparity and gt- error
				distances.push_back(dist2plane * TO_MM);
				disparities.push_back(bf_factor / point.length() - bf_factor / plane_intersect.length());
				// The negative dist2plane represents a point closer to the camera than the fitted plane
				if (ground_truth_mm) gt_errors.push_back(plane_fit_to_ground_truth_mm + (dist2plane * TO_MM));
			}

			// Show Z accuracy metric only when Ground Truth is available
			//z_accuracy->enable(ground_truth_mm > 0);
			if (ground_truth_mm > 0)
			{
				std::sort(begin(gt_errors), end(gt_errors));
				auto gt_median = gt_errors[gt_errors.size() / 2];
				auto accuracy = TO_PERCENT * (gt_median / ground_truth_mm);
				//z_accuracy->add_value(accuracy);
				//if (record) samples.push_back({ z_accuracy->get_name(),  accuracy });
				DepthQuality.ZAccuracy = accuracy;
			}
			else
			{
				DepthQuality.ZAccuracy = -1;
			}

			// Calculate Sub-pixel RMS for Stereo-based Depth sensors
			double total_sq_disparity_diff = 0;
			for (auto disparity : disparities)
			{
				total_sq_disparity_diff += disparity*disparity;
			}
			auto rms_subpixel_val = static_cast<float>(std::sqrt(total_sq_disparity_diff / disparities.size()));
			//sub_pixel_rms_error->add_value(rms_subpixel_val);
			//if (record) samples.push_back({ sub_pixel_rms_error->get_name(),  rms_subpixel_val });

			// Calculate Plane Fit RMS  (Spatial Noise) mm
			double plane_fit_err_sqr_sum = std::inner_product(distances.begin(), distances.end(), distances.begin(), 0.);
			auto rms_error_val = static_cast<float>(std::sqrt(plane_fit_err_sqr_sum / distances.size()));
			auto rms_error_val_per = TO_PERCENT * (rms_error_val / distance_mm);
			//plane_fit_rms_error->add_value(rms_error_val_per);
			//if (record) samples.push_back({ plane_fit_rms_error->get_name(),  rms_error_val });

			DepthQuality.PlaneFitError = rms_error_val_per;
			DepthQuality.SubpixelError = rms_subpixel_val;
		}
	);

	DepthQuality.bPlaneFit = rs2::is_valid(metrics.plane_corners);
	DepthQuality.Distance = metrics.distance;
	DepthQuality.Angle = metrics.angle;

	if (DepthQuality.bPlaneFit)
	{
		const float sx = 1.0f / (float)depth_profile.width();
		const float sy = 1.0f / (float)depth_profile.height();

		DepthQuality.RoiMin = FVector2D(metrics.roi.min_x * sx, metrics.roi.min_y * sy);
		DepthQuality.RoiMax = FVector2D(metrics.roi.max_x * sx, metrics.roi.max_y * sy);
	}
}

void ARealSenseInspector::ResetDepthQuality()
{
	DepthQuality.bPlaneFit = false;
	DepthQuality.Distance = -1;
	DepthQuality.Angle = -1;
	DepthQuality.FillRate = -1;
	DepthQuality.PlaneFitError = -1;
	DepthQuality.SubpixelError = -1;
	DepthQuality.ZAccuracy = -1;
	DepthQuality.RoiMin = FVector2D(0, 0);
	DepthQuality.RoiMax = FVector2D(0, 0);
}

void ARealSenseInspector::UpdatePointCloud()
{
	SCOPED_PROFILER;

	const size_t NumPoints = RsPoints->size();
	const size_t DensityPoints = (size_t)FMath::RoundToInt(NumPoints * FMath::Clamp(PclDensity, 0.0f, 1.0f));
	if (!DensityPoints) { PclMesh->SetVisibility(false); return;  }

	const size_t Step = (size_t)FMath::RoundToInt(NumPoints / (float)DensityPoints);
	if (!Step) { PclMesh->SetVisibility(false); return; }

	const size_t RenderPoints = (NumPoints / Step);
	const size_t RenderVertices = RenderPoints * 4;
	const size_t RenderIndices = RenderPoints * 6;
	const size_t MaxBufferIndices = (MAX_BUFFER_U16 / 6) * 6;
	const size_t NumSections = (RenderIndices / MaxBufferIndices) + (RenderIndices % MaxBufferIndices ? 1 : 0);

	#if 1
	REALSENSE_TRACE(TEXT(
		"PCL Id=%zu NumPoints=%zu "
		"RenderPoints=%zu RenderVertices=%zu RenderIndices=%zu NumSections=%zu "
		"Density=%.3f Step=%zu"), 
		PclFramesetId, NumPoints,
		RenderPoints, RenderVertices, RenderIndices, NumSections,
		PclDensity, Step
	);
	#endif

	const rs2::vertex* SrcVertices = RsPoints->get_vertices();
	const rs2::texture_coordinate* SrcTexcoords = RsPoints->get_texture_coordinates();

	FPointCloudVertex PV[4];
	const float Size = (PclVoxelSize * 0.5f);
	size_t SectionId = 0;
	size_t PointId = 0;
	size_t VertexId = 0;
	size_t IndexId = 0;
	size_t NumInvalid = 0;

	while (PointId < NumPoints && SectionId < NumSections)
	{
		if (!PclMeshData.Contains(SectionId))
		{
			NAMED_PROFILER("AllocMeshSection");
			PclMeshData.Add(SectionId, MakeUnique<FMeshSection>());
			PclMeshData[SectionId]->PclVertices.SetNumUninitialized(MAX_BUFFER_U16, false);
			PclMeshData[SectionId]->PclIndices.SetNumUninitialized(MAX_BUFFER_U16, false);
		}

		if (!PclMesh->DoesSectionExist(SectionId))
		{
			NAMED_PROFILER("CreateMeshSection");
			PclMesh->CreateMeshSection(SectionId, PclMeshData[SectionId]->PclVertices, PclMeshData[SectionId]->PclIndices, false, EUpdateFrequency::Frequent);
			PclMesh->SetMaterial(SectionId, PclMaterial);
		}

		FPointCloudVertex* DstVertices = &(PclMeshData[SectionId]->PclVertices[0]);
		int32* DstIndices = &(PclMeshData[SectionId]->PclIndices[0]);
		VertexId = 0;
		IndexId = 0;

		{NAMED_PROFILER("FillMeshBuffers");
		while (PointId < NumPoints && IndexId + 6 <= MaxBufferIndices)
		{
			const auto V = SrcVertices[PointId];
			if (!V.z)
			{
				PointId += Step;
				NumInvalid++;
				continue;
			}

			// the positive x-axis points to the right
			// the positive y-axis points down
			// the positive z-axis points forward
			const FVector Pos = FVector(V.z, V.x, -V.y) * 100.0f * PclScale; // meters to mm

			PV[0].Position = FVector(Pos.X, Pos.Y - Size, Pos.Z - Size);
			PV[1].Position = FVector(Pos.X, Pos.Y - Size, Pos.Z + Size);
			PV[2].Position = FVector(Pos.X, Pos.Y + Size, Pos.Z + Size);
			PV[3].Position = FVector(Pos.X, Pos.Y + Size, Pos.Z - Size);

			PV[0].Normal = FVector(-1, 0, 0);
			PV[1].Normal = FVector(-1, 0, 0);
			PV[2].Normal = FVector(-1, 0, 0);
			PV[3].Normal = FVector(-1, 0, 0);

			PV[0].Tangent = FVector(0, 0, 1);
			PV[1].Tangent = FVector(0, 0, 1);
			PV[2].Tangent = FVector(0, 0, 1);
			PV[3].Tangent = FVector(0, 0, 1);

			const auto T = SrcTexcoords[PointId];
			PV[0].UV0 = FVector2D(T.u, T.v);
			PV[1].UV0 = FVector2D(T.u, T.v);
			PV[2].UV0 = FVector2D(T.u, T.v);
			PV[3].UV0 = FVector2D(T.u, T.v);

			DstVertices[0] = PV[0];
			DstVertices[1] = PV[1];
			DstVertices[2] = PV[2];
			DstVertices[3] = PV[3];

			DstIndices[0] = VertexId + 0;
			DstIndices[1] = VertexId + 2;
			DstIndices[2] = VertexId + 1;
			DstIndices[3] = VertexId + 0;
			DstIndices[4] = VertexId + 3;
			DstIndices[5] = VertexId + 2;

			DstVertices += 4;
			DstIndices += 6;
			VertexId += 4;
			IndexId += 6;
			PointId += Step;
		}}

		if (IndexId > 0)
		{
			#if 0
			REALSENSE_TRACE(TEXT("section id=%zu points=%zu vert=%zu ind=%zu total=%zu invalid=%zu"),
				SectionId, VertexId / 4, VertexId, IndexId, PointId, NumInvalid);
			#endif

			NAMED_PROFILER("UpdateMeshSection");
			PclMeshData[SectionId]->PclVertices.SetNumUninitialized(VertexId, false);
			PclMeshData[SectionId]->PclIndices.SetNumUninitialized(IndexId, false);
			PclMesh->UpdateMeshSection(SectionId, PclMeshData[SectionId]->PclVertices, PclMeshData[SectionId]->PclIndices); // TODO: SLOW AS HELL
			PclMesh->SetMeshSectionVisible(SectionId, true);

			SectionId++;
		}
	}

	const auto NumMeshSections = PclMesh->GetNumSections();
	if (NumMeshSections > SectionId)
	{
		for (size_t i = SectionId; i < NumMeshSections; ++i)
		{
			//REALSENSE_TRACE(TEXT("hide id=%zu"), i);
			PclMesh->SetMeshSectionVisible(i, false);
		}
	}

	//REALSENSE_TRACE(TEXT("total=%zu invalid=%zu"), PointId, NumInvalid);
	PclMesh->SetVisibility(SectionId != 0);
}

void ARealSenseInspector::SetPointCloudMaterial(int SectionId, UMaterialInterface* Material)
{
	PclMaterial = Material;
	if (PclMesh)
	{
		for (auto i = 0; i < PclMesh->GetNumSections(); ++i)
		{
			PclMesh->SetMaterial(i, Material);
		}
	}
}

void ARealSenseInspector::GetDeviceInfo(rs2::device* Device)
{
	auto dpt_sensor = std::make_shared<rs2::depth_sensor>(Device->first<rs2::depth_sensor>());

	DepthScale = dpt_sensor->get_depth_scale();
	REALSENSE_TRACE(TEXT("DepthScale=%f"), DepthScale);

	auto baseline_mm = -1.f;
	auto profiles = dpt_sensor->get_stream_profiles();
	auto right_sensor = std::find_if(profiles.begin(), profiles.end(), [](rs2::stream_profile& p)
	{ return (p.stream_index() == 2) && (p.stream_type() == RS2_STREAM_INFRARED); });

	if (right_sensor != profiles.end())
	{
		auto left_sensor = std::find_if(profiles.begin(), profiles.end(), [](rs2::stream_profile& p)
		{ return (p.stream_index() == 0) && (p.stream_type() == RS2_STREAM_DEPTH); });

		auto extrin = (*left_sensor).get_extrinsics_to(*right_sensor);
		baseline_mm = fabs(extrin.translation[0]) * 1000;  // baseline in mm
	}

	StereoBaseline = baseline_mm;
	REALSENSE_TRACE(TEXT("StereoBaseline=%f"), StereoBaseline);
}

void ARealSenseInspector::EnsureProfileSupported(URealSenseDevice* Device, ERealSenseStreamType StreamType, ERealSenseFormatType Format, FRealSenseStreamMode Mode)
{
	FRealSenseStreamProfile Profile;
	Profile.StreamType = StreamType;
	Profile.Format = Format;
	Profile.Width = Mode.Width;
	Profile.Height = Mode.Height;
	Profile.Rate = Mode.Rate;

	if (!Device->SupportsProfile(Profile))
	{
		throw std::runtime_error("Profile not supported");
	}
}
