#include "PCH.h"

FDynamicTexture::FDynamicTexture(FString Name, 
	int Width, int Height, EPixelFormat Format, TextureCompressionSettings Compression) 
{
	SCOPED_PROFILER;

	REALSENSE_TRACE(TEXT("+FDynamicTexture %p %s"), this, *Name);

	this->Name = Name;
	this->Width = Width;
	this->Height = Height;
	this->Bpp = GPixelFormats[Format].BlockBytes;
	this->Format = Format;
	this->Compression = Compression;

	CommandCounter.Increment();
	ENQUEUE_UNIQUE_RENDER_COMMAND_ONEPARAMETER(
		CreateTextureCmd,
		FDynamicTexture*, Context, this,
	{
		Context->RenderCmd_CreateTexture();
	});
}

FDynamicTexture::~FDynamicTexture()
{
	SCOPED_PROFILER;

	REALSENSE_TRACE(TEXT("~FDynamicTexture %p %s CommandCounter=%d TudPool=%d"), this, *Name, (int)CommandCounter.GetValue(), (int)TudPool.size());

	{
		NAMED_PROFILER("WaitCommandCounter");
		const float Timeout = 5;
		const float Granularity = 200.0f / 1000.0f;
		int NumTries = (int)(Timeout / Granularity);

		while (CommandCounter.GetValue() > 0 && NumTries > 0)
		{
			FPlatformProcess::Sleep(Granularity);
			NumTries--;
		}

		if (CommandCounter.GetValue())
		{
			REALSENSE_ERR(TEXT("~FDynamicTexture %p %s CommandCounter=%d"), this, *Name, (int)CommandCounter.GetValue());
		}
	}

	if (TextureObject)
	{
		TextureObject->RemoveFromRoot();
		TextureObject = nullptr;
	}

	DeallocBuffer(EnqueTud);

	{
		FScopeLock Lock(&TudMx);
		for (FTextureUpdateData* Tud : TudPool)
		{
			if (Tud->Data) FMemory::Free(Tud->Data);
			FMemory::Free(Tud);
			BufferCounter.Decrement();
		}
		TudPool.clear();
	}

	if (BufferCounter.GetValue())
	{
		REALSENSE_ERR(TEXT("~FDynamicTexture %p %s BufferCounter=%d"), this, *Name, (int)BufferCounter.GetValue());
	}
}

void FDynamicTexture::Tick_GameThread()
{
	SCOPED_PROFILER;

	FTextureUpdateData* Tud = nullptr;
	if (EnqueTud)
	{
		FScopeLock Lock(&EnqueMx);
		Tud = EnqueTud;
		EnqueTud = nullptr;
	}

	// render commands should be submitted only from game thread because UE4 is braindamaged :D
	if (Tud)
	{
		CommandCounter.Increment();
		ENQUEUE_UNIQUE_RENDER_COMMAND_ONEPARAMETER(
			UpdateTextureCmd,
			FTextureUpdateData*, Tud, Tud,
		{
			Tud->Context->RenderCmd_UpdateTexture(Tud);
		});
	}
}

//
// UPDATE
//

void FDynamicTexture::EnqueUpdateData(FTextureUpdateData* Tud)
{
	SCOPED_PROFILER;

	FTextureUpdateData* Old = nullptr;
	{
		FScopeLock Lock(&EnqueMx);
		Old = EnqueTud;
		EnqueTud = Tud;
	}
	DeallocBuffer(Old);
}

void FDynamicTexture::EnqueUpdateFrame(const rs2::video_frame& Frame)
{
	SCOPED_PROFILER;

	const auto fw = Frame.get_width();
	const auto fh = Frame.get_height();
	const auto fbpp = Frame.get_bytes_per_pixel();

	if (Width != fw || Height != fh || Bpp != fbpp)
	{
		REALSENSE_ERR(TEXT("Invalid video_frame: %s Width=%d Height=%d Bpp=%d"),
			*Name, fw, fh, fbpp);
		return;
	}

	auto* Tud = AllocBuffer();
	CopyData(Tud, Frame);
	EnqueUpdateData(Tud);
}

void FDynamicTexture::CopyData(FTextureUpdateData* Tud, const rs2::video_frame& Frame)
{
	SCOPED_PROFILER;

	const auto fw = Frame.get_width();
	const auto fh = Frame.get_height();
	const auto fbpp = Frame.get_bytes_per_pixel();
	const auto fs = Frame.get_stride_in_bytes();

	uint8_t* Dst = (uint8_t*)Tud->Data;
	const uint8_t* Src = (const uint8_t*)Frame.get_data();

	for (int y = 0; y < fh; ++y)
	{
		FMemory::Memcpy(Dst, Src, fw * fbpp);
		Dst += (fw * fbpp);
		Src += fs;
	}
}

//
// TUD POOL
//

FTextureUpdateData* FDynamicTexture::AllocBuffer()
{
	SCOPED_PROFILER;

	if (!TudPool.empty())
	{
		FScopeLock Lock(&TudMx);
		if (!TudPool.empty())
		{
			auto Tud = TudPool.back();
			TudPool.pop_back();
			return Tud;
		}
	}

	auto Tud = (FTextureUpdateData*)FMemory::Malloc(sizeof(FTextureUpdateData), PLATFORM_CACHE_LINE_SIZE);
	BufferCounter.Increment();
	Tud->Context = this;
	Tud->DataSize = Width * Height * Bpp;
	Tud->Data = FMemory::Malloc(Tud->DataSize, PLATFORM_CACHE_LINE_SIZE);
	Tud->Stride = Width * Bpp;
	Tud->Width = Width;
	Tud->Height = Height;

	return Tud;
}

void FDynamicTexture::DeallocBuffer(FTextureUpdateData* Tud)
{
	if (Tud)
	{
		FScopeLock Lock(&TudMx);
		TudPool.push_back(Tud);
	}
}

//
// RENDER COMMANDS
//

void FDynamicTexture::RenderCmd_CreateTexture()
{
	SCOPED_PROFILER;

	REALSENSE_TRACE(TEXT("RenderCmd_CreateTexture %p %s"), this, *Name);

	auto Tex = UTexture2D::CreateTransient(Width, Height, Format);
	if (!Tex)
	{
		REALSENSE_ERR(TEXT("UTexture2D::CreateTransient failed"));
	}
	else
	{
		#if WITH_EDITORONLY_DATA
		Tex->MipGenSettings = TextureMipGenSettings::TMGS_NoMipmaps;
		#endif

		Tex->CompressionSettings = Compression;
		Tex->Filter = TextureFilter::TF_Nearest;
		Tex->SRGB = 0;

		Tex->AddToRoot();
		Tex->UpdateResource();

		TextureObject = Tex;
	}

	CommandCounter.Decrement();
}

void FDynamicTexture::RenderCmd_UpdateTexture(FTextureUpdateData* Tud)
{
	SCOPED_PROFILER;

	auto Tex = Tud->Context->TextureObject;
	if (Tex && Tex->Resource)
	{
		RHIUpdateTexture2D(
			((FTexture2DResource*)Tex->Resource)->GetTexture2DRHI(),
			0,
			FUpdateTextureRegion2D(0, 0, 0, 0, Tud->Width, Tud->Height),
			Tud->Stride,
			(const uint8*)Tud->Data
		);
	}

	DeallocBuffer(Tud);

	CommandCounter.Decrement();
}
