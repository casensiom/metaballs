#ifndef _METABALLS_H_
#define _METABALLS_H_

#include "marching_cube.h"

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

// PUBLIC METHODS

Grid grid_create(Index3d count, Vector3d size, Vector3d pos);
void grid_destroy(Grid *grid);

void push_ball(BallArray *balls, Ball ball);
void move_balls(Grid grid, BallArray balls, float t);
Ball random_ball(Grid grid);

void generate_mesh(Grid *grid, BallArray balls);

#endif //  _METABALLS_H_

#ifdef METABALLS_IMPLEMENTATION

static size_t cache_hit = 0;
static size_t cache_miss = 0;

Vector3d vector3d_normalize(Vector3d v)
{
    float magnitude = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
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
    return pos.x +
           (pos.y * grid.count.x) +
           (pos.z * grid.count.x * grid.count.y);
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
    Voxel voxel;

    // voxel.triangles.count = 0;
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
        Vector3d computedPoint1 = compute_edge_point(marching_cube_triangles[voxel->bits][t + 0], world_pos, voxel, *grid);
        Vector3d computedPoint2 = compute_edge_point(marching_cube_triangles[voxel->bits][t + 1], world_pos, voxel, *grid);
        Vector3d computedPoint3 = compute_edge_point(marching_cube_triangles[voxel->bits][t + 2], world_pos, voxel, *grid);

        Triangle triangle;
        triangle.vertex1 = compute_vertex(computedPoint1, *grid, balls);
        triangle.vertex2 = compute_vertex(computedPoint2, *grid, balls);
        triangle.vertex3 = compute_vertex(computedPoint3, *grid, balls);

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
    if ((voxel.neighbors & 2) != 0 && origin.x - 1 > 0)
        queue_voxel_push((Index3d){origin.x - 1, origin.y + 0, origin.z + 0}, grid);

    if ((voxel.neighbors & 4) != 0 && origin.y + 1 < grid->count.y - 1)
        queue_voxel_push((Index3d){origin.x + 0, origin.y + 1, origin.z + 0}, grid);
    if ((voxel.neighbors & 8) != 0 && origin.y - 1 > 0)
        queue_voxel_push((Index3d){origin.x + 0, origin.y - 1, origin.z + 0}, grid);

    if ((voxel.neighbors & 16) != 0 && origin.z + 1 < grid->count.z - 1)
        queue_voxel_push((Index3d){origin.x + 0, origin.y + 0, origin.z + 1}, grid);
    if ((voxel.neighbors & 32) != 0 && origin.z - 1 > 0)
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

    if (orig.x < 0 || orig.x >= grid.count.x - 1 ||
        orig.y < 0 || orig.y >= grid.count.y - 1 ||
        orig.z < 0 || orig.z >= grid.count.z - 1)
    {
        return false;
    }

    SEARCH_FIELD_LIMIT_LOOP(pos.x >= 0 && pos.x < grid.count.x - 1, pos.x--);
    SEARCH_FIELD_LIMIT_LOOP(pos.x < grid.count.x - 1, pos.x++);
    SEARCH_FIELD_LIMIT_LOOP(pos.y >= 0 && pos.y < grid.count.y - 1, pos.y--);
    SEARCH_FIELD_LIMIT_LOOP(pos.y < grid.count.y - 1, pos.y++);
    SEARCH_FIELD_LIMIT_LOOP(pos.z >= 0 && pos.z < grid.count.z - 1, pos.z--);
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
            .x = grid.center.x + radius_x * cos(balls.items[i].angular_speed.x * t + balls.items[i].initial_pos.x),
            .y = grid.center.y + radius_y * sin(balls.items[i].angular_speed.y * t + balls.items[i].initial_pos.y),
            .z = grid.center.z + radius_z * cos(balls.items[i].angular_speed.z * t + balls.items[i].initial_pos.z)};
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
        .x = 1.0f - (float)rand() / ((float)RAND_MAX / 2.0f),
        .y = 1.0f - (float)rand() / ((float)RAND_MAX / 2.0f),
        .z = 1.0f - (float)rand() / ((float)RAND_MAX / 2.0f)};

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

    float ball_mass = 50 + rand() / (RAND_MAX / 200);

    return (Ball){.pos = p, .mass = ball_mass, .color = color, .angular_speed = speed, .initial_pos = initial};
}

#endif // METABALLS_IMPLEMENTATION
