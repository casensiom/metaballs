#include "raylib_wrapper.h"

#define METABALLS_IMPLEMENTATION
#include "../metaballs.h"

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
    DrawSphere(center, 2, BLACK);

    // for (size_t i = 0; i < (grid.mesh.vertexCount) * 3; i++)
    // {
    //     Vector3 orig = {.x = 0,
    //                     .y = 0,
    //                     .z = 0};
    //     Vector3 orig2 = {.x = grid.mesh.vertices[i + 0],
    //                     .y = grid.mesh.vertices[i + 1],
    //                     .z = grid.mesh.vertices[i + 2]};

    //     Vector3 normal = {  .x = grid.mesh.normals[i + 0],
    //                         .y = grid.mesh.normals[i + 1],
    //                         .z = grid.mesh.normals[i + 2]};
    //     Vector3 end = { .x = orig.x + normal.x * grid.size.x,
    //                     .y = orig.y + normal.y * grid.size.y,
    //                     .z = orig.z + normal.z * grid.size.z};

    //     DrawLine3D(orig, end, BLUE);
    // }

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

    //         DrawLine3D(orig, end, PINK);
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
