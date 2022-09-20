#pragma once
#include <MyMath.hpp>
#include <Scene.hpp>

class Szen : public Scene
{
public:
	void Update() override;

private:
	int numberCircle = 4;
	void DrawSzen(int iteration);
};
