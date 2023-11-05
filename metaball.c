#include "marching_cube.h"

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <raylib.h>
#include <rlgl.h>

#define THRESHOLD (0.3)

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
    double x;
    double y;
    double z;
} Vector3d;

typedef struct vector2_struct
{
    double x;
    double y;
} Vector2d;

typedef struct vertex3_struct
{
    Vector3d pos;
    Vector3d normal;
    Vector2d text_coord; /* u, v, textureid */
    Color color;
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
    double mass;
    Color color;

    // orbital params
    Vector3d angular_speed;
    Vector3d initial_pos;
} Ball;
MB_DEFINE_ARRAY(Ball);

typedef size_t Index;
MB_DEFINE_ARRAY(Index);

typedef double Energy;
MB_DEFINE_ARRAY(Energy);

typedef struct cache_struct
{
    bool dirty;
    Vector3d value;
} Cache;

typedef struct node_struct
{
    double energy;
    Cache cache[MB_COORD_COUNT];
    bool visited;
} Node;
MB_DEFINE_ARRAY(Node);

typedef struct voxelcorner_struct
{
    Index3d grid_pos;
    size_t grid_index;
    double energy;
} VoxelCorner;

typedef struct triangle_struct
{
    Vertex3d vertex1;
    Vertex3d vertex2;
    Vertex3d vertex3;
} Triangle;
MB_DEFINE_FIXED_ARRAY(Triangle, 5);

typedef struct voxel_struct
{
    VoxelCorner corners[8];
    size_t bits;
    uint8_t neighbors;
    TriangleArray5 triangles;
} Voxel;

typedef struct grid_struct
{
    Vector3d pos;
    Index3d count; // number of subdivision of space
    Vector3d size; // size of separation between space subdivisions
    Vector3d center;
    NodeArray nodes;
    Index3dArray queue; // list of pending voxels to compute
} Grid;

static size_t cache_hit = 0;
static size_t cache_miss = 0;
static Mesh *metaball_mesh = 0x0;

Vector3d
vector3d_normalize(Vector3d v)
{
    double magnitude = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    v.x /= magnitude;
    v.y /= magnitude;
    v.z /= magnitude;
    return v;
}

Vector3d
vector3d_cross_product(Vector3d v1, Vector3d v2)
{
    Vector3d result;
    result.x = v1.y * v2.z - v1.z * v2.y;
    result.y = v1.z * v2.x - v1.x * v2.z;
    result.z = v1.x * v2.y - v1.y * v2.x;
    return result;
}

double
vector3d_dot_product(Vector3d v1, Vector3d v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vector3d
compute_grid_point_coordinates(Index3d pos, Grid grid)
{
    Vector3d out;
    out.x = (pos.x * grid.size.x);
    out.y = (pos.y * grid.size.y);
    out.z = (pos.z * grid.size.z);
    return out;
}

size_t
compute_grid_point_index(Index3d pos, Grid grid)
{
    return pos.x +
           (pos.y * grid.count.x) +
           (pos.z * grid.count.x * grid.count.y);
}

Mesh init_mesh(Grid grid)
{
    Mesh mesh = {0};

    size_t max_triangles = grid.count.x * grid.count.y * grid.count.z * 5; // maximum 5 triangles per voxel
    size_t max_vertex = max_triangles * 3;

    mesh.vertexCount = max_vertex;
    mesh.triangleCount = max_triangles;
    mesh.vertices = (float *)RL_MALLOC(max_vertex * 3 * sizeof(float));
    mesh.texcoords = (float *)RL_MALLOC(max_vertex * 2 * sizeof(float));
    mesh.colors = (char *)RL_MALLOC(max_vertex * 4 * sizeof(char));
    mesh.normals = (float *)RL_MALLOC(max_vertex * 3 * sizeof(float));
    mesh.indices = (unsigned short *)RL_MALLOC(max_vertex * sizeof(unsigned short));

    UploadMesh(&mesh, true);

    mesh.vertexCount = 0;
    mesh.triangleCount = 0;

    return mesh;
}

void draw_mesh(Mesh *mesh, Material material, Material wired, Matrix transform)
{

    UpdateMeshBuffer(*mesh, 0, mesh->vertices, mesh->vertexCount * 3 * sizeof(float), 0);
    UpdateMeshBuffer(*mesh, 1, mesh->texcoords, mesh->vertexCount * 2 * sizeof(float), 0);
    UpdateMeshBuffer(*mesh, 2, mesh->normals, mesh->vertexCount * 3 * sizeof(float), 0);
    UpdateMeshBuffer(*mesh, 3, mesh->colors, mesh->vertexCount * 4 * sizeof(unsigned char), 0);
    UpdateMeshBuffer(*mesh, 6, mesh->indices, mesh->triangleCount * 3 * sizeof(unsigned short), 0);

    DrawMesh(*mesh, material, transform);


    rlEnableWireMode();
    DrawMesh(*mesh, wired, transform);
    rlDisableWireMode();

    mesh->vertexCount = 0;
    mesh->triangleCount = 0;
}

void release_mesh(Mesh *mesh)
{
    UnloadMesh(*mesh);

    mesh->vertices = 0x0;
    mesh->texcoords = 0x0;
    mesh->colors = 0x0;
    mesh->normals = 0x0;
    mesh->indices = 0x0;

    mesh->vertexCount = 0;
    mesh->triangleCount = 0;
}

void mesh_add_vertex(Mesh *mesh, Vertex3d v)
{
    // Mesh vertices position array
    mesh->vertices[3 * mesh->vertexCount + 0] = v.pos.x;
    mesh->vertices[3 * mesh->vertexCount + 1] = v.pos.y;
    mesh->vertices[3 * mesh->vertexCount + 2] = v.pos.z;

    // Mesh texcoords array
    mesh->texcoords[2 * mesh->vertexCount + 0] = v.text_coord.x;
    mesh->texcoords[2 * mesh->vertexCount + 1] = v.text_coord.y;

    // Mesh colors array
    mesh->colors[4 * mesh->vertexCount + 0] = v.color.r;
    mesh->colors[4 * mesh->vertexCount + 1] = v.color.g;
    mesh->colors[4 * mesh->vertexCount + 2] = v.color.b;
    mesh->colors[4 * mesh->vertexCount + 3] = v.color.a;

    // Mesh normals array
    mesh->normals[3 * mesh->vertexCount + 0] = v.normal.x;
    mesh->normals[3 * mesh->vertexCount + 1] = v.normal.y;
    mesh->normals[3 * mesh->vertexCount + 2] = v.normal.z;

    mesh->indices[mesh->vertexCount] = mesh->vertexCount;
    mesh->vertexCount++;
}

void mesh_add_triangle(Mesh *mesh, Triangle t)
{
    // Mesh indices array initialization
    mesh_add_vertex(mesh, t.vertex1);
    mesh_add_vertex(mesh, t.vertex2);
    mesh_add_vertex(mesh, t.vertex3);

    mesh->triangleCount++;
}

void set_energy_point(Index3d grid_pos, Grid grid, double energy)
{
    size_t grid_index = compute_grid_point_index(grid_pos, grid);
    grid.nodes.items[grid_index].energy = energy;
    grid.nodes.items[grid_index].cache[MB_COORD_X].dirty = true;
    grid.nodes.items[grid_index].cache[MB_COORD_Y].dirty = true;
    grid.nodes.items[grid_index].cache[MB_COORD_Z].dirty = true;
    grid.nodes.items[grid_index].visited = false;
}

double
compute_energy(Vector3d world_pos, BallArray balls)
{
    //  e = mass/distance^2

    size_t i;
    double energy = 0;
    for (i = 0; i < balls.count; i++)
    {
        double dist_x = world_pos.x - balls.items[i].pos.x;
        double dist_y = world_pos.y - balls.items[i].pos.y;
        double dist_z = world_pos.z - balls.items[i].pos.z;
        double dist = (dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);

        if (dist < 0.0001f)
        {
            dist = 0.0001f;
        }

        energy += balls.items[i].mass / dist;
        // energy += 1.0f / sqrtf(dist);
    }

    return energy;
}

Vertex3d
compute_vertex(Vector3d vertex, Grid grid, BallArray balls)
{
    Vertex3d v;
    v.pos = vertex;

    // Compute vertex normal
    // n += 2 * mass * vector / distance^4
    Vector3d normal = {0};
    double total_energy = 0;
    for (size_t i = 0; i < balls.count; i++)
    {
        Vector3d dir = {
            .x = v.pos.x - balls.items[i].pos.x,
            .y = v.pos.y - balls.items[i].pos.y,
            .z = v.pos.z - balls.items[i].pos.z};

        double dist = dir.x * dir.x + dir.y * dir.y + dir.z * dir.z;
        total_energy += balls.items[i].mass / dist;

        normal.x += 2 * balls.items[i].mass * dir.x / (dist * dist);
        normal.y += 2 * balls.items[i].mass * dir.y / (dist * dist);
        normal.z += 2 * balls.items[i].mass * dir.z / (dist * dist);
    }

    Color color = {0};
    for (size_t i = 0; i < balls.count; i++)
    {
        Vector3d dir = {
            .x = v.pos.x - balls.items[i].pos.x,
            .y = v.pos.y - balls.items[i].pos.y,
            .z = v.pos.z - balls.items[i].pos.z};

        double dist = dir.x * dir.x + dir.y * dir.y + dir.z * dir.z;
        double energy = balls.items[i].mass / dist;
        double factor = energy / total_energy;

        color.r += (balls.items[i].color.r * factor);
        color.g += (balls.items[i].color.g * factor);
        color.b += (balls.items[i].color.b * factor);
    }
    color.a = 0xff;

    // normalize vector
    normal = vector3d_normalize(normal);

    // Calculate the spherical coordinates
    double theta = atan2(normal.z, normal.x);
    double phi = acos(normal.y);

    // Convert spherical coordinates to texture coordinates
    v.text_coord.x = (theta + PI) / (2.0 * PI);
    v.text_coord.y = phi / PI;

    v.color = color;

    return v;
}

Voxel compute_voxel(Index3d grid_pos, Grid grid)
{
    Voxel voxel;

    voxel.triangles.count = 0;
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

        voxel.bits |= voxel.corners[c].energy > THRESHOLD ? (1 << c) : 0;
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
        double t = (THRESHOLD - voxel->corners[index0].energy) / (voxel->corners[index1].energy - voxel->corners[index0].energy);

        double field_value[MB_COORD_COUNT];
        for (size_t i = 0; i < MB_COORD_COUNT; ++i)
        {
            double field_intensity_0 = marching_cube_vertices[index0][i];
            double field_intensity_1 = marching_cube_vertices[index1][i];
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

void compute_voxel_triangles(Voxel *voxel, Grid grid, BallArray balls)
{

    // This table indicates what axis variations have the second index of marching_cube_edges
    size_t coords[12] = {MB_COORD_X, MB_COORD_Z, MB_COORD_X, MB_COORD_Z,
                         MB_COORD_X, MB_COORD_Z, MB_COORD_X, MB_COORD_Z,
                         MB_COORD_Y, MB_COORD_Y, MB_COORD_Y, MB_COORD_Y};

    Index3d grid_pos = voxel->corners[0].grid_pos;
    Vector3d world_pos = compute_grid_point_coordinates(grid_pos, grid);

    Vector3d vertices[3];
    size_t i = 0;

    for (size_t t = 0; marching_cube_triangles[voxel->bits][t] != -1; t += 3)
    {
        Vector3d computedPoint1 = compute_edge_point(marching_cube_triangles[voxel->bits][t + 0], world_pos, voxel, grid);
        Vector3d computedPoint2 = compute_edge_point(marching_cube_triangles[voxel->bits][t + 1], world_pos, voxel, grid);
        Vector3d computedPoint3 = compute_edge_point(marching_cube_triangles[voxel->bits][t + 2], world_pos, voxel, grid);

        Triangle triangle;
        triangle.vertex1 = compute_vertex(computedPoint1, grid, balls);
        triangle.vertex2 = compute_vertex(computedPoint2, grid, balls);
        triangle.vertex3 = compute_vertex(computedPoint3, grid, balls);

        voxel->triangles.items[voxel->triangles.count] = triangle;
        voxel->triangles.count++;
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
                double energy = compute_energy(world_pos, balls);
                set_energy_point(grid_pos, grid, energy);
            }
        }
    }
}

void build_triangle_mesh(Voxel voxel)
{
    for (size_t i = 0; i < voxel.triangles.count; i++)
    {
        mesh_add_triangle(metaball_mesh, voxel.triangles.items[i]);
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

Index3d
queue_voxel_pop(Grid *grid)
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

void generate_mesh(Camera camera, Grid *grid, BallArray balls)
{

    for (size_t i = 0; i < balls.count; i++)
    {
        Vector3d pos = balls.items[i].pos;
        Index3d pos_voxel;
        pos_voxel.x = (size_t)floor(pos.x / grid->size.x);
        pos_voxel.y = (size_t)floor(pos.y / grid->size.y);
        pos_voxel.z = (size_t)floor(pos.z / grid->size.z);

        if (!search_field_limit(pos_voxel, &pos_voxel, *grid))
        {
            continue;
        }

        queue_voxel_push(pos_voxel, grid);
        while (!queue_voxel_empty(grid))
        {
            pos_voxel = queue_voxel_pop(grid);

            Voxel voxel = compute_voxel(pos_voxel, *grid);
            compute_voxel_triangles(&voxel, *grid, balls);
            build_triangle_mesh(voxel);
            queue_voxel_push_neighbors(voxel, grid);
        }
    }
}

void render(Grid grid)
{
    Vector3 center = (Vector3){
        .x = grid.center.x,
        .y = grid.center.y,
        .z = grid.center.z};
    Vector3 size = (Vector3){
        .x = grid.count.x * grid.size.x,
        .y = grid.count.y * grid.size.y,
        .z = grid.count.z * grid.size.z};
    DrawCubeWiresV(center, size, BLACK);
}

Grid make_grid(Index3d count, Vector3d size, Vector3d pos)
{
    Grid grid;

    double center_x = pos.x + ((count.x * size.x) + 0.5) / 2.0;
    double center_y = pos.y + ((count.y * size.y) + 0.5) / 2.0;
    double center_z = pos.z + ((count.z * size.z) + 0.5) / 2.0;

    Vector3d center = {
        .x = center_x,
        .y = center_y,
        .z = center_z};

    grid.count = count;
    grid.size = size;
    grid.pos = pos;
    grid.center = center;
    grid.nodes = MB_CREATE_ARRAY(Node, grid.count.x * grid.count.y * grid.count.z);
    grid.queue = MB_CREATE_ARRAY(Index3d, grid.count.x * grid.count.y * grid.count.z);

    return grid;
}

void free_grid(Grid grid)
{
    grid.count = (Index3d){.x = 0, .y = 0, .z = 0};
    grid.size = (Vector3d){.x = 0, .y = 0, .z = 0};
    grid.pos = (Vector3d){.x = 0, .y = 0, .z = 0};
    grid.center = (Vector3d){.x = 0, .y = 0, .z = 0};

    MB_DESTROY_ARRAY(grid.nodes);
    MB_DESTROY_ARRAY(grid.queue);
}

void push_ball(BallArray *balls, Ball ball)
{
    balls->items[balls->count] = ball;
    balls->count++;
}

void move_balls(Grid grid, BallArray balls, double t)
{
    double radius_x = (grid.count.x * grid.size.x) / 4;
    double radius_y = (grid.count.y * grid.size.y) / 4;
    double radius_z = (grid.count.z * grid.size.z) / 4;

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
    double radius_x = (grid.count.x * grid.size.x) / 4;
    double radius_y = (grid.count.y * grid.size.y) / 4;
    double radius_z = (grid.count.z * grid.size.z) / 4;

    Vector3d speed = {
        .x = 0.1 + (float)rand() / ((float)RAND_MAX / 2.0),
        .y = 0.1 + (float)rand() / ((float)RAND_MAX / 2.0),
        .z = 0.1 + (float)rand() / ((float)RAND_MAX / 2.0)};

    Vector3d initial = {
        .x = 1.0 - (float)rand() / ((float)RAND_MAX / 2.0),
        .y = 1.0 - (float)rand() / ((float)RAND_MAX / 2.0),
        .z = 1.0 - (float)rand() / ((float)RAND_MAX / 2.0)};

    Vector3d p =  (Vector3d){
        .x = grid.center.x + radius_x * cos(initial.x),
        .y = grid.center.y + radius_y * sin(initial.y),
        .z = grid.center.z + radius_z * cos(initial.z)};

    Color pallete[] = {LIGHTGRAY, GRAY, DARKGRAY, YELLOW, GOLD, ORANGE, PINK, RED,
                       MAROON, GREEN, LIME, DARKGREEN, SKYBLUE, BLUE, DARKBLUE, PURPLE,
                       VIOLET, DARKPURPLE, BEIGE, BROWN, DARKBROWN};

    size_t color_index = rand() / (RAND_MAX / 20);
    double ball_mass = 50 + rand() / (RAND_MAX / 200);

    return (Ball){.pos = p, .mass = ball_mass, .color = pallete[color_index], .initial_pos = initial, .angular_speed = speed};
}

int main(int argc, char *argv[])
{

    Index3d count = {.x = 20, .y = 20, .z = 20};
    Vector3d size = {.x = 10, .y = 10, .z = 10};
    Vector3d pos = {.x = 0, .y = 0, .z = 0};
    Grid grid = make_grid(count, size, pos);

    BallArray balls = MB_CREATE_ARRAY(Ball, 5);
    for (size_t i = 0; i < balls.capacity; i++)
    {
        push_ball(&balls, random_ball(grid));
    }

    const int screenWidth = 640;
    const int screenHeight = 480;

    // Define the camera to look into our 3d world
    Camera camera = {0};
    camera.position = (Vector3){grid.center.x, grid.center.y, grid.center.z - grid.count.z * grid.size.z * 2}; // Camera position
    camera.target = (Vector3){grid.center.x, grid.center.y, grid.center.z};                                    // Camera looking at point
    camera.up = (Vector3){0.0f, 1.0f, 0.0f};                                                                   // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                                                                       // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;                                                                    // Camera projection type

    InitWindow(screenWidth, screenHeight, "Metaballs");


    Material wired = LoadMaterialDefault();
    wired.maps[MATERIAL_MAP_DIFFUSE].color = GRAY;

    Material material = LoadMaterialDefault();
    material.maps[MATERIAL_MAP_DIFFUSE].color = WHITE;
    // material.maps[MATERIAL_MAP_DIFFUSE].texture = LoadTexture("../texel_checker.png");

    Matrix transform = (Matrix){1.0f, 0.0f, 0.0f, 0.0f,
                                0.0f, 1.0f, 0.0f, 0.0f,
                                0.0f, 0.0f, 1.0f, 0.0f,
                                0.0f, 0.0f, 0.0f, 1.0f};

    Mesh mesh = init_mesh(grid);
    metaball_mesh = &mesh;

    DisableCursor();
    while (!WindowShouldClose()) // Detect window close button or ESC key
    {
        cache_hit = 0;
        cache_miss = 0;
        // manage_keys(grid, &camera);
        UpdateCamera(&camera, CAMERA_THIRD_PERSON);
        camera.target = (Vector3){grid.center.x, grid.center.y, grid.center.z};

        if (!IsKeyDown(KEY_SPACE))
        {
            move_balls(grid, balls, GetTime());
        }

        compute_energy_field(balls, grid);
        generate_mesh(camera, &grid, balls);

        BeginDrawing();
            ClearBackground(RAYWHITE);
            BeginMode3D(camera);
                draw_mesh(metaball_mesh, material, wired, transform);
                render(grid);
            EndMode3D();

            char label_fps[256];
            sprintf(label_fps, "FPS: %d", GetFPS());
            DrawText(label_fps, 10, 10, 20, BLACK);

            sprintf(label_fps, "CACHE HIT: %lu, MISS: %lu", cache_hit, cache_miss);
            DrawText(label_fps, 10, 50, 20, BLACK);
        EndDrawing();
    }

    release_mesh(metaball_mesh);

    CloseWindow();

    MB_DESTROY_ARRAY(balls);

    free_grid(grid);

    return 0;
}
