#include "ball.h"
#include "raylib_wrapper.h"

int main(int argc, char *argv[])
{
    const int screenWidth = 800;
    const int screenHeight = 400;
    InitWindow(screenWidth, screenHeight, "Metaballs");
    //DisableCursor();

    Index3d count = {.x = 10, .y = 10, .z = 10};
    Vector3d size = {.x = 2, .y = 2, .z = 2};
    Vector3d pos = {.x = 0, .y = 0, .z = 0};
    Grid grid = grid_create(count, size, pos);

    BallArray balls = MB_CREATE_ARRAY(Ball, 3);

    Vector3d zero = {
        .x = 0.0f,
        .y = 0.0f,
        .z = 0.0f};

    Vector3d center = (Vector3d){
        .x = grid.center.x,
        .y = grid.center.y,
        .z = grid.center.z};

    Color4C color = (Color4C){
        .r = 0xff,
        .g = 0xff,
        .b = 0xff,
        .a = 0xff
    };
    ball_push(&balls, (Ball){.pos = center, .mass = 15, .color = color});

    // Define the camera to look into our 3d world
    Camera camera = {0};
    camera.position = (Vector3){grid.center.x, grid.center.y, grid.center.z - grid.count.z * grid.size.z * 1.5}; // Camera position
    camera.target = (Vector3){grid.center.x, grid.center.y, grid.center.z};                                    // Camera looking at point
    camera.up = (Vector3){0.0f, 1.0f, 0.0f};                                                                   // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                                                                       // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;                                                                    // Camera projection type


    Material wired = LoadMaterialDefault();
    wired.maps[MATERIAL_MAP_DIFFUSE].color = GRAY;
    wired.maps[MATERIAL_MAP_SPECULAR].color = GRAY;

    Texture2D text2d = LoadTexture("assets/uv.png");
    SetTextureWrap(text2d, TEXTURE_WRAP_CLAMP);
    Material material = LoadMaterialDefault();
    material.maps[MATERIAL_MAP_DIFFUSE].texture = text2d;
    material.maps[MATERIAL_MAP_SPECULAR].color = GRAY;

    Matrix transform = (Matrix){1.0f, 0.0f, 0.0f, 0.0f,
                                0.0f, 1.0f, 0.0f, 0.0f,
                                0.0f, 0.0f, 1.0f, 0.0f,
                                0.0f, 0.0f, 0.0f, 1.0f};

    // RL upload mesh
    Mesh rayMesh = raylib_upload_mesh(grid);

    while (!WindowShouldClose()) // Detect window close button or ESC key
    {

        UpdateCamera(&camera, CAMERA_THIRD_PERSON);
        camera.target = (Vector3){grid.center.x, grid.center.y, grid.center.z};

        if(IsKeyPressed(KEY_U)) {
            grid.config.wired = !grid.config.wired;
        }
        if(IsKeyPressed(KEY_P)) {
            grid.config.mapping = MAPPING_PLANAR;
        }
        if(IsKeyPressed(KEY_L)) {
            grid.config.mapping = MAPPING_SPHERIC;
        }

        if(IsKeyPressed(KEY_H)) {
            grid.config.threshold += 0.05f;
        }
        if(IsKeyPressed(KEY_J)) {
            grid.config.threshold -= 0.05f;
        }

        generate_mesh(&grid, balls);

        BeginDrawing();
            ClearBackground(RAYWHITE);
            BeginMode3D(camera);
                // RL render mesh
                raylib_render_mesh(grid, &rayMesh, material, wired, transform);
            EndMode3D();

            char label_fps[256];
            sprintf(label_fps, "FPS: %d", GetFPS());
            DrawText(label_fps, 10, 10, 20, BLACK);

            // sprintf(label_fps, "CACHE HIT: %lu, MISS: %lu", cache_hit, cache_miss);
            // DrawText(label_fps, 10, 50, 20, BLACK);
        EndDrawing();
    }


    // RL unload mesh
    raylib_unload_mesh(&rayMesh);
    MB_DESTROY_ARRAY(balls);
    grid_destroy(&grid);

    CloseWindow();

    return 0;
}
