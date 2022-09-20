#include "raylib.h"
#include "Scene.hpp"
#include "MyMath.hpp"
#include "Interface.hpp"
#include "PseudoShader.hpp"

unsigned int UIElements::currentIndex = 0;

int main(void)
{
    const int screenWidth = 1440;
    const int screenHeight = 900;

    InitWindow(screenWidth, screenHeight, "raylib [core] example - basic window");

    SetTargetFPS(60);              
    
    Scene* scenes[1] = { new PseudoShader() };
    for (size_t i = 0; i < 1; i++)
    {
        scenes[i]->screenHeight = screenHeight;
        scenes[i]->screenWidth = screenWidth;
    }
    int sceneIndex = 0;

    float x, y, z;

    while (!WindowShouldClose())
    {
        BeginDrawing();

        ClearBackground(BLACK);
        scenes[sceneIndex]->Update();
        UIElements::Update();
        UIElements::Slider(x, 0, 1, "x");
        UIElements::Slider(y, 0, 1, "y");
        UIElements::Slider(z, 0, 1, "z");
        

        EndDrawing();

    }

    CloseWindow();       
    return 0;
}
