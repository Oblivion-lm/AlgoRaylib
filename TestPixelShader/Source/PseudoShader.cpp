#include "raylib.h"
#include "MyMath.hpp"
#include "PseudoShader.hpp"

void PseudoShader::Update()
{
	for (size_t i = 0; i < screenWidth; i++)
	{
		for (size_t j = 0; j < screenHeight; j++)
		{
			PixelShader(i, j);
		}
	}
}

void PseudoShader::PixelShader(int x, int y)
{
	DrawPixel(x, y, { (unsigned char)x, (unsigned char)y, 2, 255 });
}