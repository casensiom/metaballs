
#define METABALLS_IMPLEMENTATION
#include "metaballs.h"


#include <raylib.h>
#include <rlgl.h>


Mesh raylib_upload_mesh(Grid grid)
{
    Mesh mesh = {0};

    mesh.vertexCount = grid.mesh.vertexCapacity;
    mesh.triangleCount = grid.mesh.triangleCapacity;
    mesh.vertices = grid.mesh.vertices;
    mesh.texcoords = grid.mesh.texcoords;
    mesh.colors = grid.mesh.colors;
    mesh.normals = grid.mesh.normals;
    mesh.indices = grid.mesh.indices;

    UploadMesh(&mesh, true);

    return mesh;
}

void raylib_render_mesh(Grid grid, Mesh *mesh, Material material, Material wired, Matrix transform)
{
    mesh->vertexCount = grid.mesh.vertexCount;
    mesh->triangleCount = grid.mesh.triangleCount;

    UpdateMeshBuffer(*mesh, 0, mesh->vertices, mesh->vertexCount * 3 * sizeof(float), 0);
    UpdateMeshBuffer(*mesh, 1, mesh->texcoords, mesh->vertexCount * 2 * sizeof(float), 0);
    UpdateMeshBuffer(*mesh, 2, mesh->normals, mesh->vertexCount * 3 * sizeof(float), 0);
    UpdateMeshBuffer(*mesh, 3, mesh->colors, mesh->vertexCount * 4 * sizeof(unsigned char), 0);
    UpdateMeshBuffer(*mesh, 6, mesh->indices, mesh->triangleCount * 3 * sizeof(unsigned short), 0);

    DrawMesh(*mesh, material, transform);

    if(grid.config.wired) {
        rlEnableWireMode();
        DrawMesh(*mesh, wired, transform);
        rlDisableWireMode();
    }

    Vector3 center = (Vector3){
        .x = grid.center.x - grid.size.x * 0.5f,
        .y = grid.center.y - grid.size.y * 0.5f,
        .z = grid.center.z - grid.size.z * 0.5f};
    Vector3 size = (Vector3){
        .x = (grid.count.x-1) * grid.size.x,
        .y = (grid.count.y-1) * grid.size.y,
        .z = (grid.count.z-1) * grid.size.z};
    DrawCubeWiresV(center, size, BLACK);

    Vector3 orig = (Vector3){
        .x = grid.pos.x,
        .y = grid.pos.y,
        .z = grid.pos.z};
    DrawSphere(orig, 2, BLACK);
    
    
    // for (size_t x = 0; x < grid.count.x; x++)
    // {
    //     for (size_t z = 0; z < grid.count.z; z++)
    //     {
    //         Vector3 orig = (Vector3){
    //             .x = grid.pos.x + x * grid.size.x,
    //             .y = grid.pos.y,
    //             .z = grid.pos.z + z * grid.size.z};
    //         Vector3 end = (Vector3){
    //             .x = grid.pos.x + x * grid.size.x,
    //             .y = grid.pos.y + (grid.count.y - 1) * grid.size.y,
    //             .z = grid.pos.z + z * grid.size.z};

    //         DrawLine3D(orig, end, BLUE);
    //     }
    // }
}

void raylib_unload_mesh(Mesh *mesh)
{
    // Unload rlgl mesh vboId data
    rlUnloadVertexArray(mesh->vaoId);
    if (mesh->vboId != NULL) 
    {
        for (int i = 0; i < 7; i++) 
        {
            rlUnloadVertexBuffer(mesh->vboId[i]);
        }
    }
    free(mesh->vboId);


    mesh->vertices = 0x0;
    mesh->texcoords = 0x0;
    mesh->colors = 0x0;
    mesh->normals = 0x0;
    mesh->indices = 0x0;

    mesh->vertexCount = 0;
    mesh->triangleCount = 0;
}

int main(int argc, char *argv[])
{
    const int screenWidth = 640;
    const int screenHeight = 480;
    InitWindow(screenWidth, screenHeight, "Metaballs");
    DisableCursor();

    Index3d count = {.x = 20, .y = 20, .z = 20};
    Vector3d size = {.x = 10, .y = 10, .z = 10};
    Vector3d pos = {.x = 0, .y = 0, .z = 0};
    Grid grid = grid_create(count, size, pos);

    BallArray balls = MB_CREATE_ARRAY(Ball, 5);
    for (size_t i = 0; i < balls.capacity; i++)
    {
        push_ball(&balls, random_ball(grid));
    }

    // Define the camera to look into our 3d world
    Camera camera = {0};
    camera.position = (Vector3){grid.center.x, grid.center.y, grid.center.z - grid.count.z * grid.size.z * 2}; // Camera position
    camera.target = (Vector3){grid.center.x, grid.center.y, grid.center.z};                                    // Camera looking at point
    camera.up = (Vector3){0.0f, 1.0f, 0.0f};                                                                   // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                                                                       // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;                                                                    // Camera projection type


    Material wired = LoadMaterialDefault();
    wired.maps[MATERIAL_MAP_DIFFUSE].color = GRAY;

    Texture2D text2d = LoadTexture("assets/texel_checker.png");
    SetTextureWrap(text2d, TEXTURE_WRAP_MIRROR_REPEAT);
    Material material = LoadMaterialDefault();
    material.maps[MATERIAL_MAP_DIFFUSE].texture = text2d;

    Matrix transform = (Matrix){1.0f, 0.0f, 0.0f, 0.0f,
                                0.0f, 1.0f, 0.0f, 0.0f,
                                0.0f, 0.0f, 1.0f, 0.0f,
                                0.0f, 0.0f, 0.0f, 1.0f};

    // RL upload mesh
    Mesh rayMesh = raylib_upload_mesh(grid);

    while (!WindowShouldClose()) // Detect window close button or ESC key
    {
        cache_hit = 0;
        cache_miss = 0;

        UpdateCamera(&camera, CAMERA_THIRD_PERSON);
        camera.target = (Vector3){grid.center.x, grid.center.y, grid.center.z};

        if (!IsKeyDown(KEY_SPACE))
        {
            move_balls(grid, balls, GetTime());
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
            sprintf(label_fps, "FPS: %d, THRESHOLD: %0.2f", GetFPS(), grid.config.threshold);
            DrawText(label_fps, 10, 10, 20, BLACK);

            sprintf(label_fps, "CACHE HIT: %lu, MISS: %lu", cache_hit, cache_miss);
            DrawText(label_fps, 10, 50, 20, BLACK);
        EndDrawing();
    }


    // RL unload mesh
    raylib_unload_mesh(&rayMesh);
    MB_DESTROY_ARRAY(balls);
    grid_destroy(&grid);

    CloseWindow();

    return 0;
}
