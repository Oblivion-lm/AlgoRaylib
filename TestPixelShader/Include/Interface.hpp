#pragma once
#include <MyMath.hpp>
#include <string>
#include <raylib.h>
#include <vector>

class UIElements
{

public:

	static inline void Slider(float& value, const float min = 0, const float max = 1, const std::string& text = "")
	{
		DrawText((text + " : ").c_str(), 20, 30 + currentIndex * 40, 20, {255, 255, 255, 255});

		Vector2 mouse = GetMousePosition();
		if (IsMouseButtonDown(0) && mouse.x >= 120 && mouse.x <= 240 && mouse.y > 30 + currentIndex * 40 && mouse.y < 55 + currentIndex * 40)
		{
			value = (mouse.x - 120) / (120);
		}
		float fillingWidth = (value - min) / (max - min);
		DrawRectangle(120, 30 + currentIndex * 40, fillingWidth * 120, 25, { 0, 0, 255, 255 });

		DrawRectangleLines(120, 30 + currentIndex * 40, 120, 25, {255, 255, 255, 255});
		DrawText(TextFormat("%.2f", value), 255, 30 + currentIndex * 40, 20, { 255, 255, 255, 255 });
		currentIndex++;
	}

	static inline void Update()
	{
		currentIndex = 0;
	}

private:
	static unsigned int currentIndex;
};