#pragma once

class Scene
{
public:
	int screenWidth, screenHeight;
	virtual void Update() = 0;
};
