#pragma once

struct FTextureUpdateData
{
	class FDynamicTexture* Context;
	void* Data;
	uint32 DataSize;
	uint32 Stride;
	uint32 Width;
	uint32 Height;
};

class FDynamicTexture
{
private:

	FString Name;
	UTexture2D* TextureObject = nullptr;
	int Width = 0;
	int Height = 0;
	int Bpp = 0;
	EPixelFormat Format = PF_Unknown;
	TextureCompressionSettings Compression = TextureCompressionSettings::TC_VectorDisplacementmap;

	std::vector<FTextureUpdateData*> TudPool;
	FTextureUpdateData* EnqueTud = nullptr;

	FCriticalSection TudMx;
	FCriticalSection EnqueMx;
	FThreadSafeCounter CommandCounter;
	FThreadSafeCounter BufferCounter;

private:

	void RenderCmd_CreateTexture();
	void RenderCmd_UpdateTexture(FTextureUpdateData* Tud);

public:

	FDynamicTexture(FString Name, int Width, int Height, EPixelFormat Format, TextureCompressionSettings Compression);
	virtual ~FDynamicTexture();

	void Tick_GameThread();

	void EnqueUpdateData(FTextureUpdateData* Tud);
	void EnqueUpdateFrame(const rs2::video_frame& Frame);
	void CopyData(FTextureUpdateData* Tud, const rs2::video_frame& Frame);

	FTextureUpdateData* AllocBuffer();
	void DeallocBuffer(FTextureUpdateData* Tud);

	inline UTexture2D* GetTextureObject() { return TextureObject; }
	int GetWidth() const { return Width; }
	int GetHeight() const { return Height; }
	int GetBpp() const { return Bpp; }
	EPixelFormat GetFormat() const { return Format; }
	TextureCompressionSettings GetCompression() const { return Compression; }
};
