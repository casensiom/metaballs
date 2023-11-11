#include "ball.h"
#include "raylib_wrapper.h"


void update_balls(BallArray balls, Vector3d pos, float radius) {

    float t = GetTime();

    Vector3d movement[] = {{.x = 1.1, .y = 1.3, .z = 1.7},
                            {.x = 0.7, .y = 1.9, .z = 2.3},
                            {.x = 0.3, .y = 2.9, .z = 1.1},
                            {.x = 1.3, .y = 1.7, .z = 0.7},
                            {.x = 2.3, .y = 1.9, .z = 2.9}};
    const size_t count = sizeof(movement) / sizeof(movement[0]);

    
    float ec = radius;
    for (size_t i = 0; i < balls.count; i++)
    {
        Vector3d m = movement[i % count];
	    balls.items[i].pos = (Vector3d){.x = pos.x - cos(t*m.x) * radius, 
                                        .y = pos.y - cos(t*m.y) * radius, 
                                        .z = pos.z - cos(t*m.z) * radius};
        
    }
}


int main(int argc, char *argv[])
{
    const int screenWidth = 800;
    const int screenHeight = 400;
    InitWindow(screenWidth, screenHeight, "Metaballs");
    DisableCursor();

    Index3d count = {.x = 20 * 2, .y = 20 * 2, .z = 20 * 2};
    Vector3d size = {.x = 10 / 2, .y = 10 / 2, .z = 10 / 2};
    Vector3d pos = {.x = 0, .y = 0, .z = 0};
    Grid grid = grid_create(count, size, pos);

    BallArray balls = MB_CREATE_ARRAY(Ball, 5);
    for (size_t i = 0; i < balls.capacity; i++)
    {
        ball_push(&balls, ball_create_random(grid));
    }

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
    SetTextureWrap(text2d, TEXTURE_WRAP_MIRROR_REPEAT);
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

        if (!IsKeyDown(KEY_SPACE))
        {
            update_balls(balls, grid.center, grid.size.x * 10);
        }
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
