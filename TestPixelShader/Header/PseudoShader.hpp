#pragma once
#include "MyMath.hpp"
#include "Scene.hpp"

class PseudoShader : public Scene
{
public:
	void Update() override;

private:
	void PixelShader(int x, int y);
};
