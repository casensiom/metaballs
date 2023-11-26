#include "ball.h"
#include "raylib_wrapper.h"


void update_balls(BallArray balls, size_t w, size_t d, Grid grid) 
{
    size_t index = 0;
    for (size_t x = 0; x < w; x++)
    {
        float u = (float)x / (float)w;
        for (size_t z = 0; z < d; z++)
        {
            float wave = 10.0 + 12.0 * sin(32 * u + GetTime()  * 4.0f);
            balls.items[index++].pos.y = wave;
        }
    }

    float z = grid.pos.z + grid.size.z * (grid.count.z/2);
    float wave = 12.0 * sin(GetTime() * 2.0f);
    balls.items[balls.count-1].pos.z = z + wave;
    balls.items[balls.count-2].pos.z = z - wave;

}

int main(int argc, char *argv[])
{
    const int screenWidth = 800;
    const int screenHeight = 400;
    InitWindow(screenWidth, screenHeight, "Metaballs");
    // DisableCursor();

    Index3d count = {.x = 40, .y = 10, .z = 20};
    Vector3d size = {.x = 10, .y = 10, .z = 10};
    Vector3d pos = {.x = 0, .y = 0, .z = 0};
    Grid grid = metaball_grid_create(count, size, pos);
    grid.config.wired = true;
    grid.config.threshold = 0.2;

    size_t w = 16;
    size_t d = 4;
    BallArray balls = MB_CREATE_ARRAY(Ball, w*d + 2);

    for (size_t x = 0; x < w; x++)
    {
        for (size_t z = 0; z < d; z++)
        {
            ball_push(&balls, (Ball){
                .pos = {.x = pos.x + size.x * (x * count.x/(w - 1)), 
                        .y = pos.y + 2 * size.y, 
                        .z = pos.z + size.z * (z * count.z/(d - 1))}, 
                .mass = 48, 
                .color = {.r = 0x80 , .g = 0x80, .b = (128) + ((64+64)*x)/w , .a = 0xff}});
        }
    }

    Vector3d ball_pos = {.x = pos.x + 10 * size.x, .y = pos.y + 4 * size.y, .z = pos.z + 10 * size.z};
    ball_push(&balls, (Ball){.pos = ball_pos, .mass = -96, .color = {.r = 0xff, .g = 0xff, .b = 0xff, .a = 0xff}});

    ball_pos = (Vector3d){.x = pos.x + 30 * size.x, .y = pos.y + 4 * size.y, .z = pos.z + 10 * size.z};
    ball_push(&balls, (Ball){.pos = ball_pos, .mass = -96, .color = {.r = 0xff, .g = 0xff, .b = 0xff, .a = 0xff}});


    // Define the camera to look into our 3d world
    Camera camera = {0};
    camera.position = (Vector3){grid.center.x, grid.center.y, grid.center.z - grid.count.z * grid.size.z}; // Camera position
    camera.target = (Vector3){grid.center.x, grid.center.y, grid.center.z};                                    // Camera looking at point
    camera.up = (Vector3){0.0f, 1.0f, 0.0f};                                                                   // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                                                                       // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;                                                                    // Camera projection type


    Material wired = LoadMaterialDefault();
    wired.maps[MATERIAL_MAP_DIFFUSE].color = GRAY;
    wired.maps[MATERIAL_MAP_SPECULAR].color = GRAY;

    Material material = LoadMaterialDefault();
    material.maps[MATERIAL_MAP_DIFFUSE].color = WHITE;
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
            update_balls(balls, w, d, grid);
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
        metaball_generate_mesh(&grid, balls);

        BeginDrawing();
            ClearBackground(RAYWHITE);
            BeginMode3D(camera);
                // RL render mesh
                raylib_render_mesh(grid, &rayMesh, material, wired, transform);
            EndMode3D();
        EndDrawing();
    }


    // RL unload mesh
    raylib_unload_mesh(&rayMesh);
    MB_DESTROY_ARRAY(balls);
    metaball_grid_destroy(&grid);

    CloseWindow();

    return 0;
}
