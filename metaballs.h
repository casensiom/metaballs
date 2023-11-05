// Metaballs library generator - v1.0 - public domain
// https://github.com/casensiom/metaballs
//


#ifndef _METABALLS_H_
#define _METABALLS_H_

// #include "marching_cube.h"

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#ifndef PI
#define PI 3.14159265
#endif

#define MB_DEFINE_ARRAY(TYPE)          \
    typedef struct TYPE##_struct_array \
    {                                  \
        TYPE *items;                   \
        size_t count;                  \
        size_t capacity;               \
    } TYPE##Array;

#define MB_DEFINE_FIXED_ARRAY(TYPE, SIZE)     \
    typedef struct TYPE##_struct_array_##SIZE \
    {                                         \
        TYPE items[SIZE];                     \
        size_t count;                         \
    } TYPE##Array##SIZE;

#define MB_CREATE_ARRAY(TYPE, CAPACITY)                   \
    (TYPE##Array)                                         \
    {                                                     \
        .items = (TYPE *)malloc(CAPACITY * sizeof(TYPE)), \
        .count = 0,                                       \
        .capacity = CAPACITY                              \
    }

#define MB_DESTROY_ARRAY(arr) \
    {                         \
        free(arr.items);      \
        arr.items = NULL;     \
        arr.count = 0;        \
        arr.capacity = 0;     \
    }

enum MB_COORDINATE
{
    MB_COORD_X = 0,
    MB_COORD_Y = 1,
    MB_COORD_Z = 2,
    MB_COORD_COUNT
};

typedef struct vector3_struct
{
    float x;
    float y;
    float z;
} Vector3d;

typedef struct vector2_struct
{
    float x;
    float y;
} Vector2d;

typedef struct color_struct
{
    unsigned char r;
    unsigned char g;
    unsigned char b;
    unsigned char a;
} Color4C;

typedef struct vertex3_struct
{
    Vector3d pos;
    Vector3d normal;
    Vector2d text_coord;
    Color4C color;
} Vertex3d;

typedef struct index3_struct
{
    size_t x;
    size_t y;
    size_t z;
} Index3d;
MB_DEFINE_ARRAY(Index3d);

typedef struct ball_struct
{
    Vector3d pos;
    float mass;
    Color4C color;

    // orbital params
    Vector3d angular_speed;
    Vector3d initial_pos;
} Ball;
MB_DEFINE_ARRAY(Ball);

typedef size_t Index;
MB_DEFINE_ARRAY(Index);

typedef struct cache_struct
{
    bool dirty;
    Vector3d value;
} Cache;

typedef struct node_struct
{
    float energy;
    Cache cache[MB_COORD_COUNT];
    bool visited;
} Node;
MB_DEFINE_ARRAY(Node);

typedef struct voxelcorner_struct
{
    Index3d grid_pos;
    size_t grid_index;
    float energy;
} VoxelCorner;

typedef struct triangle_struct
{
    Vertex3d vertex1;
    Vertex3d vertex2;
    Vertex3d vertex3;
} Triangle;

typedef struct mesh3d_struct
{
    size_t vertexCapacity;
    size_t triangleCapacity;
    size_t vertexCount;
    size_t triangleCount;
    float *vertices;
    float *texcoords;
    unsigned char *colors;
    float *normals;
    unsigned short *indices;
} Mesh3d;

typedef struct voxel_struct
{
    VoxelCorner corners[8];
    size_t bits;
    uint8_t neighbors;
} Voxel;

typedef enum
{
    MAPPING_PLANAR = 0,
    MAPPING_SPHERIC
} Mapping;

typedef struct config_struct
{
    bool wired;
    Mapping mapping;
    float threshold;
} Config;

typedef struct grid_struct
{
    Vector3d pos;
    Index3d count; // number of subdivision of space
    Vector3d size; // size of separation between space subdivisions
    Vector3d center;
    NodeArray nodes;
    Index3dArray queue; // list of pending voxels to compute
    Mesh3d mesh;
    Config config;
} Grid;

#if defined(__cplusplus)
extern "C" {
#endif

Grid grid_create(Index3d count, Vector3d size, Vector3d pos);

void grid_destroy(Grid *grid);

void push_ball(BallArray *balls, Ball ball);

void move_balls(Grid grid, BallArray balls, float t);

Ball random_ball(Grid grid);

void generate_mesh(Grid *grid, BallArray balls);

#if defined(__cplusplus)
}
#endif

#endif //  _METABALLS_H_

#ifdef METABALLS_IMPLEMENTATION


static char marching_cube_triangles[256][16];
static char marching_cube_edges[12][2];
static int marching_cube_vertices[8][3];
static char marching_cube_neighbors[256];

static size_t cache_hit = 0;
static size_t cache_miss = 0;

float vector3d_length(Vector3d start, Vector3d end) {
    float x = end.x - start.x;
    float y = end.y - start.y;
    float z = end.z - start.z;
    return sqrtf(x * x + y * y + z * z);
}

Vector3d vector3d_normalize(Vector3d v)
{
    float magnitude = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    v.x /= magnitude;
    v.y /= magnitude;
    v.z /= magnitude;
    return v;
}

Vector3d vector3d_cross_product(Vector3d v1, Vector3d v2)
{
    Vector3d result;
    result.x = v1.y * v2.z - v1.z * v2.y;
    result.y = v1.z * v2.x - v1.x * v2.z;
    result.z = v1.x * v2.y - v1.y * v2.x;
    return result;
}

float vector3d_dot_product(Vector3d v1, Vector3d v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}


Vector3d compute_grid_point_coordinates(Index3d pos, Grid grid)
{
    Vector3d out;
    out.x = grid.pos.x + (pos.x * grid.size.x);
    out.y = grid.pos.y + (pos.y * grid.size.y);
    out.z = grid.pos.z + (pos.z * grid.size.z);
    return out;
}

size_t compute_grid_point_index(Index3d pos, Grid grid)
{    
    size_t ret = pos.x +
           (pos.y * grid.count.x) +
           (pos.z * grid.count.x * grid.count.y);
    assert(ret < grid.count.x * grid.count.y * grid.count.z);
    return ret;
}

Index3d compute_grid_index(Vector3d pos, Grid grid) {
    Index3d out;
    out.x = (size_t)floor((pos.x - grid.pos.x) / grid.size.x);
    out.y = (size_t)floor((pos.y - grid.pos.y) / grid.size.y);
    out.z = (size_t)floor((pos.z - grid.pos.z) / grid.size.z);
    return out;
}
//--

Mesh3d mesh_create(Grid grid)
{
    Mesh3d mesh = {0};

    size_t max_triangles = grid.count.x * grid.count.y * grid.count.z * 5; // maximum 5 triangles per voxel
    size_t max_vertex = max_triangles * 3;

    mesh.vertexCapacity = max_vertex;
    mesh.triangleCapacity = max_triangles;
    mesh.vertexCount = 0;
    mesh.triangleCount = 0;
    mesh.vertices = (float *)malloc(max_vertex * 3 * sizeof(float));
    mesh.texcoords = (float *)malloc(max_vertex * 2 * sizeof(float));
    mesh.colors = (unsigned char *)malloc(max_vertex * 4 * sizeof(unsigned char));
    mesh.normals = (float *)malloc(max_vertex * 3 * sizeof(float));
    mesh.indices = (unsigned short *)malloc(max_vertex * sizeof(unsigned short));

    return mesh;
}

void mesh_destroy(Mesh3d *mesh)
{
    free(mesh->vertices);
    free(mesh->texcoords);
    free(mesh->colors);
    free(mesh->normals);
    free(mesh->indices);

    *mesh = (Mesh3d){0};
}

void mesh_add_vertex(Mesh3d *mesh, Vertex3d v)
{
    mesh->vertices[3 * mesh->vertexCount + 0] = v.pos.x;
    mesh->vertices[3 * mesh->vertexCount + 1] = v.pos.y;
    mesh->vertices[3 * mesh->vertexCount + 2] = v.pos.z;

    mesh->texcoords[2 * mesh->vertexCount + 0] = v.text_coord.x;
    mesh->texcoords[2 * mesh->vertexCount + 1] = v.text_coord.y;

    mesh->colors[4 * mesh->vertexCount + 0] = v.color.r;
    mesh->colors[4 * mesh->vertexCount + 1] = v.color.g;
    mesh->colors[4 * mesh->vertexCount + 2] = v.color.b;
    mesh->colors[4 * mesh->vertexCount + 3] = v.color.a;

    mesh->normals[3 * mesh->vertexCount + 0] = v.normal.x;
    mesh->normals[3 * mesh->vertexCount + 1] = v.normal.y;
    mesh->normals[3 * mesh->vertexCount + 2] = v.normal.z;

    mesh->indices[mesh->vertexCount] = mesh->vertexCount;
    mesh->vertexCount++;
}

void mesh_add_triangle(Mesh3d *mesh, Triangle t)
{
    mesh_add_vertex(mesh, t.vertex1);
    mesh_add_vertex(mesh, t.vertex2);
    mesh_add_vertex(mesh, t.vertex3);

    mesh->triangleCount++;
}

//--

void set_energy_point(Index3d grid_pos, Grid grid, float energy)
{
    assert(grid_pos.x < grid.count.x);
    assert(grid_pos.y < grid.count.y);
    assert(grid_pos.z < grid.count.z);

    size_t grid_index = compute_grid_point_index(grid_pos, grid);
    grid.nodes.items[grid_index].energy = energy;
    grid.nodes.items[grid_index].cache[MB_COORD_X].dirty = true;
    grid.nodes.items[grid_index].cache[MB_COORD_Y].dirty = true;
    grid.nodes.items[grid_index].cache[MB_COORD_Z].dirty = true;
    grid.nodes.items[grid_index].visited = false;
}

float compute_energy(Vector3d world_pos, BallArray balls)
{
    //  e = mass/distance^2

    size_t i;
    float energy = 0;
    for (i = 0; i < balls.count; i++)
    {
        float dist_x = world_pos.x - balls.items[i].pos.x;
        float dist_y = world_pos.y - balls.items[i].pos.y;
        float dist_z = world_pos.z - balls.items[i].pos.z;
        float dist = (dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);

        if (dist < 0.0001f)
        {
            dist = 0.0001f;
        }

        energy += balls.items[i].mass / dist;
        // energy += 1.0f / sqrtf(dist);
    }

    return energy;
}

Vertex3d compute_vertex(Vector3d vertex, Grid grid, BallArray balls)
{
    Vertex3d v;
    v.pos = vertex;

    // Compute vertex normal
    // n += 2 * mass * vector / distance^4
    Vector3d normal = {0};
    float total_energy = 0;
    for (size_t i = 0; i < balls.count; i++)
    {
        Vector3d dir = {
            .x = v.pos.x - balls.items[i].pos.x,
            .y = v.pos.y - balls.items[i].pos.y,
            .z = v.pos.z - balls.items[i].pos.z};

        float dist = dir.x * dir.x + dir.y * dir.y + dir.z * dir.z;
        total_energy += balls.items[i].mass / dist;

        normal.x += 2 * balls.items[i].mass * dir.x / (dist * dist);
        normal.y += 2 * balls.items[i].mass * dir.y / (dist * dist);
        normal.z += 2 * balls.items[i].mass * dir.z / (dist * dist);
    }

    Color4C color = {0};
    for (size_t i = 0; i < balls.count; i++)
    {
        Vector3d dir = {
            .x = v.pos.x - balls.items[i].pos.x,
            .y = v.pos.y - balls.items[i].pos.y,
            .z = v.pos.z - balls.items[i].pos.z};

        float dist = dir.x * dir.x + dir.y * dir.y + dir.z * dir.z;
        float energy = balls.items[i].mass / dist;
        float factor = energy / total_energy;

        color.r += (balls.items[i].color.r * factor);
        color.g += (balls.items[i].color.g * factor);
        color.b += (balls.items[i].color.b * factor);
    }
    color.a = 0xff;

    normal = vector3d_normalize(normal);

    if (grid.config.mapping == MAPPING_PLANAR)
    {
        v.text_coord.x = (normal.x / 2.0) + 0.5;
        v.text_coord.y = -(normal.y / 2.0) + 0.5;
    }
    else
    {
        float theta = atan2(normal.z, normal.x);
        float phi = acos(normal.y);
        v.text_coord.x = (theta + PI) / (2.0 * PI) + 0.5;
        v.text_coord.y = phi / PI + 0.5;
    }

    v.color = color;

    return v;
}

Voxel compute_voxel(Index3d grid_pos, Grid grid)
{
    assert(grid_pos.x < grid.count.x - 1);
    assert(grid_pos.y < grid.count.y - 1);
    assert(grid_pos.z < grid.count.z - 1);

    Voxel voxel = {0};
    voxel.bits = 0;
    for (size_t c = 0; c < 8; ++c)
    {
        Index3d pos = {
            .x = grid_pos.x + marching_cube_vertices[c][0],
            .y = grid_pos.y + marching_cube_vertices[c][1],
            .z = grid_pos.z + marching_cube_vertices[c][2]};
        size_t index = compute_grid_point_index(pos, grid);

        voxel.corners[c].grid_pos = pos;
        voxel.corners[c].grid_index = index;
        voxel.corners[c].energy = grid.nodes.items[index].energy;

        voxel.bits |= voxel.corners[c].energy > grid.config.threshold ? (1 << c) : 0;
    }

    voxel.neighbors = marching_cube_neighbors[voxel.bits];
    return voxel;
}

Vector3d compute_edge_point(char edge, Vector3d orig, Voxel *voxel, Grid grid)
{
    Vector3d computedPoint;
    char index0 = marching_cube_edges[edge][0];
    char index1 = marching_cube_edges[edge][1];

    // This table indicates what axis variations have the second index of marching_cube_edges
    size_t coords[12] = {MB_COORD_X, MB_COORD_Z, MB_COORD_X, MB_COORD_Z,
                         MB_COORD_X, MB_COORD_Z, MB_COORD_X, MB_COORD_Z,
                         MB_COORD_Y, MB_COORD_Y, MB_COORD_Y, MB_COORD_Y};

    size_t coord = coords[edge];
    if (!grid.nodes.items[voxel->corners[index0].grid_index].cache[coord].dirty)
    {
        computedPoint = grid.nodes.items[voxel->corners[index0].grid_index].cache[coord].value;
        cache_hit++;
    }
    else
    {
        cache_miss++;
        float t = (grid.config.threshold - voxel->corners[index0].energy) / (voxel->corners[index1].energy - voxel->corners[index0].energy);

        float field_value[MB_COORD_COUNT];
        for (size_t i = 0; i < MB_COORD_COUNT; ++i)
        {
            float field_intensity_0 = marching_cube_vertices[index0][i];
            float field_intensity_1 = marching_cube_vertices[index1][i];
            field_value[i] = field_intensity_0 + (field_intensity_1 - field_intensity_0) * t;
        }

        computedPoint = (Vector3d){
            orig.x + field_value[MB_COORD_X] * grid.size.x,
            orig.y + field_value[MB_COORD_Y] * grid.size.y,
            orig.z + field_value[MB_COORD_Z] * grid.size.z};

        grid.nodes.items[voxel->corners[index0].grid_index].cache[coord].value = computedPoint;
        grid.nodes.items[voxel->corners[index0].grid_index].cache[coord].dirty = false;
    }
    return computedPoint;
}

void compute_voxel_triangles(Voxel *voxel, Grid *grid, BallArray balls)
{
    // This table indicates what axis variations have the second index of marching_cube_edges
    size_t coords[12] = {MB_COORD_X, MB_COORD_Z, MB_COORD_X, MB_COORD_Z,
                         MB_COORD_X, MB_COORD_Z, MB_COORD_X, MB_COORD_Z,
                         MB_COORD_Y, MB_COORD_Y, MB_COORD_Y, MB_COORD_Y};

    Index3d grid_pos = voxel->corners[0].grid_pos;
    Vector3d world_pos = compute_grid_point_coordinates(grid_pos, *grid);

    Vector3d vertices[3];
    size_t i = 0;

    for (size_t t = 0; marching_cube_triangles[voxel->bits][t] != -1; t += 3)
    {
        Vector3d p1 = compute_edge_point(marching_cube_triangles[voxel->bits][t + 0], world_pos, voxel, *grid);
        Vector3d p2 = compute_edge_point(marching_cube_triangles[voxel->bits][t + 1], world_pos, voxel, *grid);
        Vector3d p3 = compute_edge_point(marching_cube_triangles[voxel->bits][t + 2], world_pos, voxel, *grid);

        Triangle triangle;
        triangle.vertex1 = compute_vertex(p1, *grid, balls);
        triangle.vertex2 = compute_vertex(p2, *grid, balls);
        triangle.vertex3 = compute_vertex(p3, *grid, balls);

        mesh_add_triangle(&(grid->mesh), triangle);
    }
}

void compute_energy_field(BallArray balls, Grid grid)
{
    size_t x, y, z;

    for (x = 0; x < grid.count.x; x++)
    {
        for (y = 0; y < grid.count.y; y++)
        {
            for (z = 0; z < grid.count.z; z++)
            {
                Index3d grid_pos = {.x = x, .y = y, .z = z};
                Vector3d world_pos = compute_grid_point_coordinates(grid_pos, grid);
                float energy = compute_energy(world_pos, balls);
                set_energy_point(grid_pos, grid, energy);
            }
        }
    }
}

void queue_voxel_push(Index3d pos, Grid *grid)
{
    assert(grid->queue.count < grid->queue.capacity);
    assert(pos.x <= grid->count.x - 1);
    assert(pos.y <= grid->count.y - 1);
    assert(pos.z <= grid->count.z - 1);
    size_t index = compute_grid_point_index(pos, *grid);
    if (grid->nodes.items[index].visited == false)
    {
        grid->nodes.items[index].visited = true;
        grid->queue.items[grid->queue.count++] = pos;
    }
}

Index3d queue_voxel_pop(Grid *grid)
{
    assert(grid->queue.count >= 0);
    return grid->queue.items[--grid->queue.count];
}

bool queue_voxel_empty(Grid *grid)
{
    return (grid->queue.count == 0);
}

void queue_voxel_push_neighbors(Voxel voxel, Grid *grid)
{
    Index3d origin = voxel.corners[0].grid_pos;

    if ((voxel.neighbors & 1) != 0 && origin.x + 1 < grid->count.x - 1)
        queue_voxel_push((Index3d){origin.x + 1, origin.y + 0, origin.z + 0}, grid);
    if ((voxel.neighbors & 2) != 0 && origin.x >= 1)
        queue_voxel_push((Index3d){origin.x - 1, origin.y + 0, origin.z + 0}, grid);

    if ((voxel.neighbors & 4) != 0 && origin.y + 1 < grid->count.y - 1)
        queue_voxel_push((Index3d){origin.x + 0, origin.y + 1, origin.z + 0}, grid);
    if ((voxel.neighbors & 8) != 0 && origin.y >= 1)
        queue_voxel_push((Index3d){origin.x + 0, origin.y - 1, origin.z + 0}, grid);

    if ((voxel.neighbors & 16) != 0 && origin.z + 1 < grid->count.z - 1)
        queue_voxel_push((Index3d){origin.x + 0, origin.y + 0, origin.z + 1}, grid);
    if ((voxel.neighbors & 32) != 0 && origin.z >= 1)
        queue_voxel_push((Index3d){origin.x + 0, origin.y + 0, origin.z - 1}, grid);
}

#define SEARCH_FIELD_LIMIT_LOOP(COMPARE, UPDATE)       \
    {                                                  \
        Index3d pos = orig;                            \
        while (COMPARE)                                \
        {                                              \
            Voxel voxel = compute_voxel(pos, grid);    \
            if (voxel.bits == 0 || voxel.bits == 0xFF) \
            {                                          \
                UPDATE;                                \
            }                                          \
            else                                       \
            {                                          \
                *found = pos;                          \
                return true;                           \
            }                                          \
        }                                              \
    }

bool search_field_limit(Index3d orig, Index3d *found, Grid grid)
{
    if (orig.x >= grid.count.x - 1 ||
        orig.y >= grid.count.y - 1 ||
        orig.z >= grid.count.z - 1)
    {
        return false;
    }

    // pos component are size_t, can't be negative
    SEARCH_FIELD_LIMIT_LOOP(pos.x < grid.count.x - 1, pos.x--);
    SEARCH_FIELD_LIMIT_LOOP(pos.x < grid.count.x - 1, pos.x++);
    SEARCH_FIELD_LIMIT_LOOP(pos.y < grid.count.y - 1, pos.y--);
    SEARCH_FIELD_LIMIT_LOOP(pos.y < grid.count.y - 1, pos.y++);
    SEARCH_FIELD_LIMIT_LOOP(pos.z < grid.count.z - 1, pos.z--);
    SEARCH_FIELD_LIMIT_LOOP(pos.z < grid.count.z - 1, pos.z++);

    return false;
}

void generate_mesh(Grid *grid, BallArray balls)
{
    grid->mesh.triangleCount = 0;
    grid->mesh.vertexCount = 0;

    compute_energy_field(balls, *grid);
    for (size_t i = 0; i < balls.count; i++)
    {
        Vector3d pos = balls.items[i].pos;
        Index3d pos_voxel = compute_grid_index(pos, *grid);

        if (!search_field_limit(pos_voxel, &pos_voxel, *grid))
        {
            continue;
        }

        queue_voxel_push(pos_voxel, grid);
        while (!queue_voxel_empty(grid))
        {
            pos_voxel = queue_voxel_pop(grid);

            Voxel voxel = compute_voxel(pos_voxel, *grid);
            compute_voxel_triangles(&voxel, grid, balls);
            queue_voxel_push_neighbors(voxel, grid);
        }
    }
}

Grid grid_create(Index3d count, Vector3d size, Vector3d pos)
{
    Grid grid;

    Vector3d center = {
        .x = pos.x + (count.x * size.x) / 2.0f,
        .y = pos.y + (count.y * size.y) / 2.0f,
        .z = pos.z + (count.z * size.z) / 2.0f};

    grid.count = count;
    grid.size = size;
    grid.pos = pos;
    grid.center = center;
    grid.nodes = MB_CREATE_ARRAY(Node, grid.count.x * grid.count.y * grid.count.z);
    grid.queue = MB_CREATE_ARRAY(Index3d, grid.count.x * grid.count.y * grid.count.z);
    grid.mesh = mesh_create(grid);
    grid.config.threshold = 0.3;

    return grid;
}

void grid_destroy(Grid *grid)
{
    grid->count = (Index3d){.x = 0, .y = 0, .z = 0};
    grid->size = (Vector3d){.x = 0, .y = 0, .z = 0};
    grid->pos = (Vector3d){.x = 0, .y = 0, .z = 0};
    grid->center = (Vector3d){.x = 0, .y = 0, .z = 0};

    MB_DESTROY_ARRAY(grid->nodes);
    MB_DESTROY_ARRAY(grid->queue);
    mesh_destroy(&(grid->mesh));
}

void push_ball(BallArray *balls, Ball ball)
{
    balls->items[balls->count] = ball;
    balls->count++;
}

void move_balls(Grid grid, BallArray balls, float t)
{
    float radius_x = (grid.count.x * grid.size.x) / 4;
    float radius_y = (grid.count.y * grid.size.y) / 4;
    float radius_z = (grid.count.z * grid.size.z) / 4;

    for (size_t i = 0; i < balls.count; i++)
    {
        balls.items[i].pos = (Vector3d){
            .x = grid.center.x + radius_x * cosf(balls.items[i].angular_speed.x * t + balls.items[i].initial_pos.x),
            .y = grid.center.y + radius_y * sinf(balls.items[i].angular_speed.y * t + balls.items[i].initial_pos.y),
            .z = grid.center.z + radius_z * cosf(balls.items[i].angular_speed.z * t + balls.items[i].initial_pos.z)};
    }
}

Ball random_ball(Grid grid)
{
    float radius_x = (grid.count.x * grid.size.x) / 4;
    float radius_y = (grid.count.y * grid.size.y) / 4;
    float radius_z = (grid.count.z * grid.size.z) / 4;

    Vector3d speed = {
        .x = 0.1f + (float)rand() / ((float)RAND_MAX / 2.0f),
        .y = 0.1f + (float)rand() / ((float)RAND_MAX / 2.0f),
        .z = 0.1f + (float)rand() / ((float)RAND_MAX / 2.0f)};

    Vector3d initial = {
        .x = 2.0f - (float)rand() / ((float)RAND_MAX / 4.0f),
        .y = 2.0f - (float)rand() / ((float)RAND_MAX / 4.0f),
        .z = 2.0f - (float)rand() / ((float)RAND_MAX / 4.0f)};

    Vector3d p = (Vector3d){
        .x = grid.center.x + radius_x * cos(initial.x),
        .y = grid.center.y + radius_y * sin(initial.y),
        .z = grid.center.z + radius_z * cos(initial.z)};

    Color4C color = (Color4C){
        .r = (unsigned char)(64  + (rand() % 192)),
        .g = (unsigned char)(128 + (rand() % 128)),
        .b = (unsigned char)(128 + (rand() % 128)),
        .a = (unsigned char)(255)
    };

    float ball_mass = 50.0f + (float)(rand() % 200);

    return (Ball){.pos = p, .mass = ball_mass, .color = color, .angular_speed = speed, .initial_pos = initial};
}






//=============================================================================
//==============================================================================

//        +----------+
//       /|7        /|6
//      / |        / |
//     +----------+  |
//    4|  |       |5 |
//     |  |       |  |
//     |  |       |  |
//     |  +-------|--+
//     | / 3      | / 2
//     |/         |/
//     +----------+
//    0            1       

//        +----------+
//      7/|    6    /|
//      / |        /5|
//     +----------+  |
//     |  |  4    |  |11
//     |10|       |  |
//    8|  |       |9 |
//     |  +-------|--+
//     | /     2  | /
//     |/3        |/1
//     +----------+
//           0           

//==============================================================================
//==============================================================================

// Cube neighbors
//
// bit 0 : x + 1
// bit 1 : x - 1
// bit 2 : y + 1
// bit 3 : y - 1
// bit 4 : z + 1
// bit 5 : z - 1
static char marching_cube_neighbors[256] = 
{
     0, 42, 41, 43, 25, 59, 57, 59, 26, 58, 59, 59, 27, 59, 59, 51, 
    38, 46, 47, 47, 63, 63, 63, 63, 62, 62, 63, 63, 63, 63, 63, 55, 
    37, 47, 45, 47, 61, 63, 61, 63, 63, 63, 63, 63, 63, 63, 63, 55, 
    39, 47, 47, 15, 63, 63, 63, 31, 63, 63, 63, 31, 63, 63, 63, 23, 
    21, 63, 61, 63, 29, 63, 61, 63, 31, 63, 63, 63, 31, 63, 63, 55, 
    55, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 55, 
    53, 63, 61, 63, 61, 63, 60, 62, 63, 63, 63, 63, 63, 63, 62, 54, 
    55, 63, 63, 31, 63, 63, 62, 30, 63, 63, 63, 31, 63, 63, 62, 22, 
    22, 62, 63, 63, 31, 63, 63, 63, 30, 62, 63, 63, 31, 63, 63, 55, 
    54, 62, 63, 63, 63, 63, 63, 63, 62, 60, 63, 61, 63, 61, 63, 53, 
    55, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 55, 
    55, 63, 63, 31, 63, 63, 63, 31, 63, 61, 63, 29, 63, 61, 63, 21, 
    23, 63, 63, 63, 31, 63, 63, 63, 31, 63, 63, 63, 15, 47, 47, 39, 
    55, 63, 63, 63, 63, 63, 63, 63, 63, 61, 63, 61, 47, 45, 47, 37, 
    55, 63, 63, 63, 63, 63, 62, 62, 63, 63, 63, 63, 47, 47, 46, 38, 
    51, 59, 59, 27, 59, 59, 58, 26, 59, 57, 59, 25, 43, 41, 42, 0
};

// Cube vertices
static int marching_cube_vertices[8][3] = 
{
    {0,0,0},
    {1,0,0},
    {1,0,1},
    {0,0,1},
    {0,1,0},
    {1,1,0},
    {1,1,1},
    {0,1,1}
};

// This is the edges and the direction on them. They are designed so
// that edges of neighboring cubes are in the same direction.
static char marching_cube_edges[12][2] = 
{
    {0,1}, {1,2}, {3,2}, {0,3},
    {4,5}, {5,6}, {7,6}, {4,7},
    {0,4}, {1,5}, {3,7}, {2,6}
};

// This list gives the edges that the triangles in each case intersect.
static char marching_cube_triangles[256][16] = 
{ 
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 0  
    { 3,  0,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 1  
    { 9,  0,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 2  
    { 3,  1,  8,  1,  9,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 3  
    {11,  1,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 4  
    { 3,  0,  8, 11,  1,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 5  
    {11,  9,  2,  9,  0,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 6  
    { 3,  2,  8,  8,  2, 11,  8, 11,  9, -1, -1, -1, -1, -1, -1, -1}, // 7  
    { 2,  3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 8 
    { 2,  0, 10,  0,  8, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 9 
    { 0,  1,  9, 10,  2,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 10  
    { 2,  1, 10, 10,  1,  9, 10,  9,  8, -1, -1, -1, -1, -1, -1, -1}, // 11  
    { 1,  3, 11,  3, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 12  
    { 1,  0, 11, 11,  0,  8, 11,  8, 10, -1, -1, -1, -1, -1, -1, -1}, // 13  
    { 0,  3,  9,  9,  3, 10,  9, 10, 11, -1, -1, -1, -1, -1, -1, -1}, // 14  
    {11,  9,  8, 10, 11,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 15  
    { 8,  4,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 16  
    { 0,  4,  3,  4,  7,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 17  
    { 9,  0,  1,  7,  8,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 18  
    { 9,  4,  1,  1,  4,  7,  1,  7,  3, -1, -1, -1, -1, -1, -1, -1}, // 19  
    {11,  1,  2,  7,  8,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 20  
    { 7,  3,  4,  4,  3,  0, 11,  1,  2, -1, -1, -1, -1, -1, -1, -1}, // 21  
    {11,  9,  2,  2,  9,  0,  7,  8,  4, -1, -1, -1, -1, -1, -1, -1}, // 22  
    { 9,  2, 11,  7,  2,  9,  3,  2,  7,  4,  7,  9, -1, -1, -1, -1}, // 23  
    { 7,  8,  4,  2,  3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 24  
    { 7, 10,  4,  4, 10,  2,  4,  2,  0, -1, -1, -1, -1, -1, -1, -1}, // 25  
    { 1,  9,  0,  7,  8,  4, 10,  2,  3, -1, -1, -1, -1, -1, -1, -1}, // 26  
    {10,  4,  7, 10,  9,  4,  2,  9, 10,  1,  9,  2, -1, -1, -1, -1}, // 27  
    { 1,  3, 11, 11,  3, 10,  4,  7,  8, -1, -1, -1, -1, -1, -1, -1}, // 28  
    {11,  1, 10, 10,  1,  4,  4,  1,  0,  4,  7, 10, -1, -1, -1, -1}, // 29  
    { 8,  4,  7, 10,  9,  0, 11,  9, 10,  3, 10,  0, -1, -1, -1, -1}, // 30  
    {10,  4,  7,  9,  4, 10, 11,  9, 10, -1, -1, -1, -1, -1, -1, -1}, // 31  
    { 4,  9,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 32  
    { 4,  9,  5,  3,  0,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 33  
    { 4,  0,  5,  0,  1,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 34  
    { 4,  8,  5,  5,  8,  3,  5,  3,  1, -1, -1, -1, -1, -1, -1, -1}, // 35  
    {11,  1,  2,  4,  9,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 36  
    { 8,  3,  0, 11,  1,  2,  5,  4,  9, -1, -1, -1, -1, -1, -1, -1}, // 37  
    {11,  5,  2,  2,  5,  4,  2,  4,  0, -1, -1, -1, -1, -1, -1, -1}, // 38  
    { 5,  2, 11,  5,  3,  2,  4,  3,  5,  8,  3,  4, -1, -1, -1, -1}, // 39  
    { 4,  9,  5, 10,  2,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 40  
    { 2,  0, 10, 10,  0,  8,  5,  4,  9, -1, -1, -1, -1, -1, -1, -1}, // 41  
    { 4,  0,  5,  5,  0,  1, 10,  2,  3, -1, -1, -1, -1, -1, -1, -1}, // 42  
    { 5,  2,  1,  8,  2,  5, 10,  2,  8,  5,  4,  8, -1, -1, -1, -1}, // 43  
    {10, 11,  3,  3, 11,  1,  4,  9,  5, -1, -1, -1, -1, -1, -1, -1}, // 44  
    { 5,  4,  9,  1,  0,  8,  1,  8, 11, 11,  8, 10, -1, -1, -1, -1}, // 45  
    { 0,  5,  4, 10,  5,  0, 11,  5, 10,  3, 10,  0, -1, -1, -1, -1}, // 46  
    { 8,  5,  4, 11,  5,  8, 10, 11,  8, -1, -1, -1, -1, -1, -1, -1}, // 47  
    { 8,  9,  7,  9,  5,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 48  
    { 0,  9,  3,  3,  9,  5,  3,  5,  7, -1, -1, -1, -1, -1, -1, -1}, // 49  
    { 8,  0,  7,  7,  0,  1,  7,  1,  5, -1, -1, -1, -1, -1, -1, -1}, // 50  
    { 3,  1,  5,  7,  3,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 51  
    { 8,  9,  7,  7,  9,  5,  2, 11,  1, -1, -1, -1, -1, -1, -1, -1}, // 52  
    { 2, 11,  1,  0,  9,  5,  0,  5,  3,  3,  5,  7, -1, -1, -1, -1}, // 53  
    { 2,  8,  0,  5,  8,  2,  7,  8,  5,  2, 11,  5, -1, -1, -1, -1}, // 54  
    { 5,  2, 11,  3,  2,  5,  7,  3,  5, -1, -1, -1, -1, -1, -1, -1}, // 55  
    { 5,  7,  9,  9,  7,  8,  2,  3, 10, -1, -1, -1, -1, -1, -1, -1}, // 56  
    { 7,  9,  5,  2,  9,  7,  0,  9,  2, 10,  2,  7, -1, -1, -1, -1}, // 57  
    {10,  2,  3,  8,  0,  1,  8,  1,  7,  7,  1,  5, -1, -1, -1, -1}, // 58  
    { 1, 10,  2,  7, 10,  1,  5,  7,  1, -1, -1, -1, -1, -1, -1, -1}, // 59  
    { 8,  9,  5,  7,  8,  5,  3, 11,  1, 10, 11,  3, -1, -1, -1, -1}, // 60  
    { 0,  5,  7,  9,  5,  0,  0,  7, 10, 11,  1,  0,  0, 10, 11, -1}, // 61  
    { 0, 10, 11,  3, 10,  0,  0, 11,  5,  7,  8,  0,  0,  5,  7, -1}, // 62  
    { 5, 10, 11,  5,  7, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 63  
    { 5, 11,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 64  
    { 3,  0,  8,  6,  5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 65  
    { 1,  9,  0,  6,  5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 66  
    { 3,  1,  8,  8,  1,  9,  6,  5, 11, -1, -1, -1, -1, -1, -1, -1}, // 67  
    { 5,  1,  6,  1,  2,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 68  
    { 5,  1,  6,  6,  1,  2,  8,  3,  0, -1, -1, -1, -1, -1, -1, -1}, // 69  
    { 5,  9,  6,  6,  9,  0,  6,  0,  2, -1, -1, -1, -1, -1, -1, -1}, // 70  
    { 8,  5,  9,  2,  5,  8,  6,  5,  2,  8,  3,  2, -1, -1, -1, -1}, // 71  
    {10,  2,  3,  5, 11,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 72  
    { 8, 10,  0,  0, 10,  2,  5, 11,  6, -1, -1, -1, -1, -1, -1, -1}, // 73  
    { 9,  0,  1, 10,  2,  3,  6,  5, 11, -1, -1, -1, -1, -1, -1, -1}, // 74  
    { 6,  5, 11,  2,  1,  9,  2,  9, 10, 10,  9,  8, -1, -1, -1, -1}, // 75  
    {10,  6,  3,  3,  6,  5,  3,  5,  1, -1, -1, -1, -1, -1, -1, -1}, // 76  
    {10,  0,  8,  5,  0, 10,  1,  0,  5,  6,  5, 10, -1, -1, -1, -1}, // 77  
    { 6,  3, 10,  6,  0,  3,  5,  0,  6,  9,  0,  5, -1, -1, -1, -1}, // 78  
    { 9,  6,  5, 10,  6,  9,  8, 10,  9, -1, -1, -1, -1, -1, -1, -1}, // 79  
    { 6,  5, 11,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 80  
    { 0,  4,  3,  3,  4,  7, 11,  6,  5, -1, -1, -1, -1, -1, -1, -1}, // 81  
    { 0,  1,  9,  6,  5, 11,  7,  8,  4, -1, -1, -1, -1, -1, -1, -1}, // 82  
    { 5, 11,  6,  7,  1,  9,  3,  1,  7,  4,  7,  9, -1, -1, -1, -1}, // 83  
    { 2,  6,  1,  1,  6,  5,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1}, // 84  
    { 5,  1,  2,  6,  5,  2,  4,  3,  0,  7,  3,  4, -1, -1, -1, -1}, // 85  
    { 7,  8,  4,  5,  9,  0,  5,  0,  6,  6,  0,  2, -1, -1, -1, -1}, // 86  
    { 9,  7,  3,  4,  7,  9,  9,  3,  2,  6,  5,  9,  9,  2,  6, -1}, // 87  
    { 2,  3, 10,  4,  7,  8,  5, 11,  6, -1, -1, -1, -1, -1, -1, -1}, // 88 
    { 6,  5, 11,  2,  4,  7,  0,  4,  2, 10,  2,  7, -1, -1, -1, -1}, // 89  
    { 9,  0,  1,  8,  4,  7, 10,  2,  3,  6,  5, 11, -1, -1, -1, -1}, // 90  
    { 1,  9,  2,  2,  9, 10, 10,  9,  4,  4,  7, 10,  6,  5, 11, -1}, // 91  
    { 7,  8,  4,  5,  3, 10,  1,  3,  5,  6,  5, 10, -1, -1, -1, -1}, // 92  
    {10,  5,  1,  6,  5, 10, 10,  1,  0,  4,  7, 10, 10,  0,  4, -1}, // 93 
    { 9,  0,  5,  5,  0,  6,  6,  0,  3,  3, 10,  6,  7,  8,  4, -1}, // 94  
    { 9,  6,  5, 10,  6,  9,  9,  4,  7,  9,  7, 10, -1, -1, -1, -1}, // 95  
    { 9, 11,  4, 11,  6,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 96  
    { 6,  4, 11, 11,  4,  9,  3,  0,  8, -1, -1, -1, -1, -1, -1, -1}, // 97  
    { 1, 11,  0,  0, 11,  6,  0,  6,  4, -1, -1, -1, -1, -1, -1, -1}, // 98 
    { 1,  8,  3,  6,  8,  1,  4,  8,  6, 11,  6,  1, -1, -1, -1, -1}, // 99 
    { 9,  1,  4,  4,  1,  2,  4,  2,  6, -1, -1, -1, -1, -1, -1, -1}, // 100  
    { 8,  3,  0,  9,  1,  2,  9,  2,  4,  4,  2,  6, -1, -1, -1, -1}, // 101  
    { 4,  0,  2,  6,  4,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 102  
    { 2,  8,  3,  4,  8,  2,  6,  4,  2, -1, -1, -1, -1, -1, -1, -1}, // 103  
    { 9, 11,  4,  4, 11,  6,  3, 10,  2, -1, -1, -1, -1, -1, -1, -1}, // 104  
    { 2,  0,  8, 10,  2,  8, 11,  4,  9,  6,  4, 11, -1, -1, -1, -1}, // 105  
    { 2,  3, 10,  6,  0,  1,  4,  0,  6, 11,  6,  1, -1, -1, -1, -1}, // 106  
    { 1,  6,  4, 11,  6,  1,  1,  4,  8, 10,  2,  1,  1,  8, 10, -1}, // 107  
    { 4,  9,  6,  6,  9,  3,  3,  9,  1,  3, 10,  6, -1, -1, -1, -1}, // 108  
    { 1,  8, 10,  0,  8,  1,  1, 10,  6,  4,  9,  1,  1,  6,  4, -1}, // 109  
    { 6,  3, 10,  0,  3,  6,  4,  0,  6, -1, -1, -1, -1, -1, -1, -1}, // 110  
    { 8,  6,  4,  8, 10,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 111  
    { 6,  7, 11, 11,  7,  8, 11,  8,  9, -1, -1, -1, -1, -1, -1, -1}, // 112  
    { 3,  0,  7,  7,  0, 11, 11,  0,  9, 11,  6,  7, -1, -1, -1, -1}, // 113  
    { 7, 11,  6,  7,  1, 11,  8,  1,  7,  0,  1,  8, -1, -1, -1, -1}, // 114  
    { 7, 11,  6,  1, 11,  7,  3,  1,  7, -1, -1, -1, -1, -1, -1, -1}, // 115 
    { 6,  1,  2,  8,  1,  6,  9,  1,  8,  7,  8,  6, -1, -1, -1, -1}, // 116  
    { 9,  2,  6,  1,  2,  9,  9,  6,  7,  3,  0,  9,  9,  7,  3, -1}, // 117  
    { 0,  7,  8,  6,  7,  0,  2,  6,  0, -1, -1, -1, -1, -1, -1, -1}, // 118  
    { 2,  7,  3,  2,  6,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 119  
    {10,  2,  3,  8, 11,  6,  9, 11,  8,  7,  8,  6, -1, -1, -1, -1}, // 120  
    { 7,  2,  0, 10,  2,  7,  7,  0,  9, 11,  6,  7,  7,  9, 11, -1}, // 121 
    { 0,  1,  8,  8,  1,  7,  7,  1, 11, 11,  6,  7, 10,  2,  3, -1}, // 122  
    { 1, 10,  2,  7, 10,  1,  1, 11,  6,  1,  6,  7, -1, -1, -1, -1}, // 123  
    { 6,  8,  9,  7,  8,  6,  6,  9,  1,  3, 10,  6,  6,  1,  3, -1}, // 124  
    { 1,  0,  9,  7, 10,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 125  
    { 0,  7,  8,  6,  7,  0,  0,  3, 10,  0, 10,  6, -1, -1, -1, -1}, // 126  
    { 6,  7, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 127  
    {10,  7,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 128  
    { 8,  3,  0,  6, 10,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 129  
    { 9,  0,  1,  6, 10,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 130  
    { 9,  8,  1,  1,  8,  3,  6, 10,  7, -1, -1, -1, -1, -1, -1, -1}, // 131  
    { 2, 11,  1,  7,  6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 132  
    {11,  1,  2,  8,  3,  0,  7,  6, 10, -1, -1, -1, -1, -1, -1, -1}, // 133  
    { 0,  2,  9,  9,  2, 11,  7,  6, 10, -1, -1, -1, -1, -1, -1, -1}, // 134  
    { 7,  6, 10,  3,  2, 11,  3, 11,  8,  8, 11,  9, -1, -1, -1, -1}, // 135  
    { 3,  7,  2,  7,  6,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 136  
    { 8,  7,  0,  0,  7,  6,  0,  6,  2, -1, -1, -1, -1, -1, -1, -1}, // 137  
    { 6,  2,  7,  7,  2,  3,  9,  0,  1, -1, -1, -1, -1, -1, -1, -1}, // 138  
    { 2,  1,  6,  6,  1,  8,  8,  1,  9,  6,  8,  7, -1, -1, -1, -1}, // 139  
    { 6, 11,  7,  7, 11,  1,  7,  1,  3, -1, -1, -1, -1, -1, -1, -1}, // 140  
    { 6, 11,  7, 11,  1,  7,  7,  1,  8,  8,  1,  0, -1, -1, -1, -1}, // 141  
    { 7,  0,  3, 11,  0,  7,  9,  0, 11,  7,  6, 11, -1, -1, -1, -1}, // 142  
    {11,  7,  6,  8,  7, 11,  9,  8, 11, -1, -1, -1, -1, -1, -1, -1}, // 143  
    { 4,  6,  8,  6, 10,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 144  
    {10,  3,  6,  6,  3,  0,  6,  0,  4, -1, -1, -1, -1, -1, -1, -1}, // 145  
    {10,  8,  6,  6,  8,  4,  1,  9,  0, -1, -1, -1, -1, -1, -1, -1}, // 146  
    { 6,  9,  4,  3,  9,  6,  1,  9,  3,  6, 10,  3, -1, -1, -1, -1}, // 147  
    { 4,  6,  8,  8,  6, 10,  1,  2, 11, -1, -1, -1, -1, -1, -1, -1}, // 148  
    {11,  1,  2, 10,  3,  0, 10,  0,  6,  6,  0,  4, -1, -1, -1, -1}, // 149  
    { 8,  4, 10, 10,  4,  6,  9,  0,  2,  9,  2, 11, -1, -1, -1, -1}, // 150  
    { 3, 11,  9,  2, 11,  3,  3,  9,  4,  6, 10,  3,  3,  4,  6, -1}, // 151  
    { 3,  8,  2,  2,  8,  4,  2,  4,  6, -1, -1, -1, -1, -1, -1, -1}, // 152  
    { 2,  0,  4,  2,  4,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 153  
    { 0,  1,  9,  4,  2,  3,  6,  2,  4,  8,  4,  3, -1, -1, -1, -1}, // 154  
    { 4,  1,  9,  2,  1,  4,  6,  2,  4, -1, -1, -1, -1, -1, -1, -1}, // 155  
    { 3,  8,  1,  1,  8,  6,  6,  8,  4,  1,  6, 11, -1, -1, -1, -1}, // 156  
    { 0, 11,  1,  6, 11,  0,  4,  6,  0, -1, -1, -1, -1, -1, -1, -1}, // 157  
    { 3,  4,  6,  8,  4,  3,  3,  6, 11,  9,  0,  3,  3, 11,  9, -1}, // 158  
    { 4, 11,  9,  4,  6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 159  
    { 5,  4,  9, 10,  7,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 160  
    { 3,  0,  8,  5,  4,  9,  6, 10,  7, -1, -1, -1, -1, -1, -1, -1}, // 161  
    { 1,  5,  0,  0,  5,  4, 10,  7,  6, -1, -1, -1, -1, -1, -1, -1}, // 162  
    { 6, 10,  7,  4,  8,  3,  4,  3,  5,  5,  3,  1, -1, -1, -1, -1}, // 163  
    { 4,  9,  5,  2, 11,  1, 10,  7,  6, -1, -1, -1, -1, -1, -1, -1}, // 164  
    { 7,  6, 10, 11,  1,  2,  3,  0,  8,  5,  4,  9, -1, -1, -1, -1}, // 165  
    {10,  7,  6, 11,  5,  4, 11,  4,  2,  2,  4,  0, -1, -1, -1, -1}, // 166  
    { 8,  3,  4,  4,  3,  5,  5,  3,  2,  2, 11,  5,  6, 10,  7, -1}, // 167  
    { 3,  7,  2,  2,  7,  6,  9,  5,  4, -1, -1, -1, -1, -1, -1, -1}, // 168  
    { 4,  9,  5,  6,  0,  8,  2,  0,  6,  7,  6,  8, -1, -1, -1, -1}, // 169  
    { 2,  3,  6,  6,  3,  7,  0,  1,  5,  0,  5,  4, -1, -1, -1, -1}, // 170  
    { 8,  6,  2,  7,  6,  8,  8,  2,  1,  5,  4,  8,  8,  1,  5, -1}, // 171 
    { 4,  9,  5,  6, 11,  1,  6,  1,  7,  7,  1,  3, -1, -1, -1, -1}, // 172  
    {11,  1,  6,  6,  1,  7,  7,  1,  0,  0,  8,  7,  4,  9,  5, -1}, // 173  
    {11,  4,  0,  5,  4, 11, 11,  0,  3,  7,  6, 11, 11,  3,  7, -1}, // 174  
    {11,  7,  6,  8,  7, 11, 11,  5,  4, 11,  4,  8, -1, -1, -1, -1}, // 175  
    { 5,  6,  9,  9,  6, 10,  9, 10,  8, -1, -1, -1, -1, -1, -1, -1}, // 176  
    {10,  3,  6,  3,  0,  6,  6,  0,  5,  5,  0,  9, -1, -1, -1, -1}, // 177  
    { 8,  0, 10, 10,  0,  5,  5,  0,  1, 10,  5,  6, -1, -1, -1, -1}, // 178  
    { 3,  6, 10,  5,  6,  3,  1,  5,  3, -1, -1, -1, -1, -1, -1, -1}, // 179  
    {11,  1,  2, 10,  9,  5,  8,  9, 10,  6, 10,  5, -1, -1, -1, -1}, // 180  
    { 3,  0, 10, 10,  0,  6,  6,  0,  9,  9,  5,  6, 11,  1,  2, -1}, // 181  
    { 5, 10,  8,  6, 10,  5,  5,  8,  0,  2, 11,  5,  5,  0,  2, -1}, // 182  
    { 3,  6, 10,  5,  6,  3,  3,  2, 11,  3, 11,  5, -1, -1, -1, -1}, // 183  
    { 9,  5,  8,  8,  5,  2,  2,  5,  6,  2,  3,  8, -1, -1, -1, -1}, // 184  
    { 6,  9,  5,  0,  9,  6,  2,  0,  6, -1, -1, -1, -1, -1, -1, -1}, // 185  
    { 8,  1,  5,  0,  1,  8,  8,  5,  6,  2,  3,  8,  8,  6,  2, -1}, // 186  
    { 6,  1,  5,  6,  2,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 187  
    { 6,  1,  3, 11,  1,  6,  6,  3,  8,  9,  5,  6,  6,  8,  9, -1}, // 188  
    { 0, 11,  1,  6, 11,  0,  0,  9,  5,  0,  5,  6, -1, -1, -1, -1}, // 189  
    { 8,  0,  3, 11,  5,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 190  
    { 6, 11,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 191  
    {11, 10,  5, 10,  7,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 192  
    {11, 10,  5,  5, 10,  7,  0,  8,  3, -1, -1, -1, -1, -1, -1, -1}, // 193  
    { 7,  5, 10, 10,  5, 11,  0,  1,  9, -1, -1, -1, -1, -1, -1, -1}, // 194  
    { 5, 11,  7,  7, 11, 10,  1,  9,  8,  1,  8,  3, -1, -1, -1, -1}, // 195  
    { 2, 10,  1,  1, 10,  7,  1,  7,  5, -1, -1, -1, -1, -1, -1, -1}, // 196  
    { 3,  0,  8,  7,  1,  2,  5,  1,  7, 10,  7,  2, -1, -1, -1, -1}, // 197  
    { 5,  9,  7,  7,  9,  2,  2,  9,  0,  7,  2, 10, -1, -1, -1, -1}, // 198  
    { 2,  7,  5, 10,  7,  2,  2,  5,  9,  8,  3,  2,  2,  9,  8, -1}, // 199  
    {11,  2,  5,  5,  2,  3,  5,  3,  7, -1, -1, -1, -1, -1, -1, -1}, // 200  
    { 0,  8,  2,  2,  8,  5,  5,  8,  7,  5, 11,  2, -1, -1, -1, -1}, // 201  
    { 1,  9,  0,  3,  5, 11,  7,  5,  3,  2,  3, 11, -1, -1, -1, -1}, // 202  
    { 2,  9,  8,  1,  9,  2,  2,  8,  7,  5, 11,  2,  2,  7,  5, -1}, // 203  
    { 5,  1,  3,  5,  3,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 204  
    { 7,  0,  8,  1,  0,  7,  5,  1,  7, -1, -1, -1, -1, -1, -1, -1}, // 205  
    { 3,  9,  0,  5,  9,  3,  7,  5,  3, -1, -1, -1, -1, -1, -1, -1}, // 206  
    { 7,  9,  8,  7,  5,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 207  
    { 4,  5,  8,  8,  5, 11,  8, 11, 10, -1, -1, -1, -1, -1, -1, -1}, // 208  
    { 4,  5,  0,  0,  5, 10, 10,  5, 11,  0, 10,  3, -1, -1, -1, -1}, // 209  
    { 9,  0,  1, 11,  8,  4, 10,  8, 11,  5, 11,  4, -1, -1, -1, -1}, // 210  
    { 4, 11, 10,  5, 11,  4,  4, 10,  3,  1,  9,  4,  4,  3,  1, -1}, // 211  
    { 1,  2,  5,  5,  2,  8,  8,  2, 10,  8,  4,  5, -1, -1, -1, -1}, // 212  
    {10,  0,  4,  3,  0, 10, 10,  4,  5,  1,  2, 10, 10,  5,  1, -1}, // 213  
    { 5,  0,  2,  9,  0,  5,  5,  2, 10,  8,  4,  5,  5, 10,  8, -1}, // 214  
    { 5,  9,  4,  3,  2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 215  
    {11,  2,  5,  2,  3,  5,  5,  3,  4,  4,  3,  8, -1, -1, -1, -1}, // 216  
    { 2,  5, 11,  4,  5,  2,  0,  4,  2, -1, -1, -1, -1, -1, -1, -1}, // 217  
    { 2,  3, 11, 11,  3,  5,  5,  3,  8,  8,  4,  5,  9,  0,  1, -1}, // 218  
    { 2,  5, 11,  4,  5,  2,  2,  1,  9,  2,  9,  4, -1, -1, -1, -1}, // 219  
    { 5,  8,  4,  3,  8,  5,  1,  3,  5, -1, -1, -1, -1, -1, -1, -1}, // 220  
    { 5,  0,  4,  5,  1,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 221  
    { 5,  8,  4,  3,  8,  5,  5,  9,  0,  5,  0,  3, -1, -1, -1, -1}, // 222  
    { 5,  9,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 223  
    { 7,  4, 10, 10,  4,  9, 10,  9, 11, -1, -1, -1, -1, -1, -1, -1}, // 224  
    { 3,  0,  8,  7,  4,  9,  7,  9, 10, 10,  9, 11, -1, -1, -1, -1}, // 225  
    {10,  1, 11,  4,  1, 10,  0,  1,  4, 10,  7,  4, -1, -1, -1, -1}, // 226  
    { 4,  3,  1,  8,  3,  4,  4,  1, 11, 10,  7,  4,  4, 11, 10, -1}, // 227  
    { 7,  4, 10,  4,  9, 10, 10,  9,  2,  2,  9,  1, -1, -1, -1, -1}, // 228  
    { 4,  9,  7,  7,  9, 10, 10,  9,  1,  1,  2, 10,  3,  0,  8, -1}, // 229  
    { 4, 10,  7,  2, 10,  4,  0,  2,  4, -1, -1, -1, -1, -1, -1, -1}, // 230  
    { 4, 10,  7,  2, 10,  4,  4,  8,  3,  4,  3,  2, -1, -1, -1, -1}, // 231  
    {11,  2,  9,  9,  2,  7,  7,  2,  3,  9,  7,  4, -1, -1, -1, -1}, // 232  
    { 7,  9, 11,  4,  9,  7,  7, 11,  2,  0,  8,  7,  7,  2,  0, -1}, // 233  
    {11,  3,  7,  2,  3, 11, 11,  7,  4,  0,  1, 11, 11,  4,  0, -1}, // 234  
    { 2,  1, 11,  4,  8,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 235  
    { 1,  4,  9,  7,  4,  1,  3,  7,  1, -1, -1, -1, -1, -1, -1, -1}, // 236  
    { 1,  4,  9,  7,  4,  1,  1,  0,  8,  1,  8,  7, -1, -1, -1, -1}, // 237  
    { 3,  4,  0,  3,  7,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 238  
    { 7,  4,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 239  
    { 8,  9, 11,  8, 11, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 240  
    { 9,  3,  0, 10,  3,  9, 11, 10,  9, -1, -1, -1, -1, -1, -1, -1}, // 241  
    {11,  0,  1,  8,  0, 11, 10,  8, 11, -1, -1, -1, -1, -1, -1, -1}, // 242  
    {11,  3,  1, 11, 10,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 243  
    {10,  1,  2,  9,  1, 10,  8,  9, 10, -1, -1, -1, -1, -1, -1, -1}, // 244  
    { 9,  3,  0, 10,  3,  9,  9,  1,  2,  9,  2, 10, -1, -1, -1, -1}, // 245  
    {10,  0,  2, 10,  8,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 246  
    {10,  3,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 247  
    { 8,  2,  3, 11,  2,  8,  9, 11,  8, -1, -1, -1, -1, -1, -1, -1}, // 248  
    { 2,  9, 11,  2,  0,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 249  
    { 8,  2,  3, 11,  2,  8,  8,  0,  1,  8,  1, 11, -1, -1, -1, -1}, // 250  
    { 2,  1, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 251  
    { 8,  1,  3,  8,  9,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 252  
    { 1,  0,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 253  
    { 8,  0,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, // 254  
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}  // 255  
};

#endif // METABALLS_IMPLEMENTATION
