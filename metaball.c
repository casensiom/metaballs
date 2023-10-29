#include "marching_cube.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <errno.h>

#include <raylib.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <profileapi.h>
#endif


#define THRESHOLD (0.3)

#define MB_DEFINE_ARRAY(TYPE) \
typedef struct TYPE ## _struct_array { \
    TYPE *items; \
    size_t count; \
    size_t capacity; \
} TYPE ## Array;

#define MB_DEFINE_FIXED_ARRAY(TYPE, SIZE) \
typedef struct TYPE ## _struct_array_ ## SIZE { \
    TYPE items[SIZE]; \
    size_t count; \
} TYPE ## Array ## SIZE;

#define MB_CREATE_ARRAY(TYPE, CAPACITY) \
(TYPE ## Array){ \
    .items = (TYPE *)malloc(CAPACITY * sizeof(TYPE)),\
    .count=0, \
    .capacity=CAPACITY \
}

#define MB_DESTROY_ARRAY(arr) \
{ \
    free(arr.items); \
    arr.items = NULL; \
    arr.count=0; \
    arr.capacity=0; \
}

enum MB_COORDINATE {
    MB_COORD_X = 0,
    MB_COORD_Y = 1,
    MB_COORD_Z = 2,
    MB_COORD_COUNT
};

typedef struct vector3_struct {
    double x;
    double y;
    double z;
} Vector3d;

typedef struct vector2_struct {
    double x;
    double y;
} Vector2d;

typedef struct vertex3_struct {
    Vector3d pos;
    Vector3d normal;
    Vector2d text_coord;    /* u, v, textureid */
} Vertex3d;

typedef struct index3_struct {
    size_t x;
    size_t y;
    size_t z;
} Index3d;

typedef struct ball_struct {
    Vector3d pos;
    double mass;
} Ball;
MB_DEFINE_ARRAY(Ball);

typedef size_t Index;
MB_DEFINE_ARRAY(Index);

typedef double Energy;
MB_DEFINE_ARRAY(Energy);

typedef struct cache_struct {
    bool dirty;
    Vector3d value;
} Cache;

typedef struct node_struct {
    double energy;
    Cache cache[MB_COORD_COUNT];
} Node;
MB_DEFINE_ARRAY(Node);

typedef struct voxelcorner_struct {
    Index3d grid_pos;
    size_t grid_index;
    double energy;
} VoxelCorner;

typedef struct triangle_struct {
    Vertex3d vertex1;
    Vertex3d vertex2;
    Vertex3d vertex3;
} Triangle;
MB_DEFINE_FIXED_ARRAY(Triangle, 5);

typedef struct voxel_struct {
    VoxelCorner corners[8];
    size_t bits;
    TriangleArray5 triangles;
} Voxel;

typedef struct grid_struct {
    Vector3d pos;
    Index3d count;  // number of subdivision of space
    Vector3d size;   // size of separation between space subdivisions
    Vector3d offset;
    NodeArray nodes;
} Grid;


Vector3d
vector3d_normalize(Vector3d v) {
    double magnitude = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    v.x /= magnitude;
    v.y /= magnitude;
    v.z /= magnitude;
    return v;
}

Vector3d
vector3d_cross_product(Vector3d v1, Vector3d v2) {
    Vector3d result;
    result.x = v1.y * v2.z - v1.z * v2.y;
    result.y = v1.z * v2.x - v1.x * v2.z;
    result.z = v1.x * v2.y - v1.y * v2.x;
    return result;
}

double
vector3d_dot_product(Vector3d v1, Vector3d v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}



Vector3d
compute_grid_point_coordinates(Index3d pos, Grid grid) {
    Vector3d out;
    out.x = (pos.x * grid.size.x) + grid.offset.x;
    out.y = (pos.y * grid.size.y) + grid.offset.y;
    out.z = (pos.z * grid.size.z) + grid.offset.z;
    return out;
}

size_t
compute_grid_point_index(Index3d pos, Grid grid) {
    return pos.x +
            (pos.y * grid.count.x) +
            (pos.z * grid.count.x * grid.count.y);
}

void
set_energy_point(Index3d grid_pos, Grid grid, double energy) {
    size_t grid_index = compute_grid_point_index(grid_pos, grid);
    grid.nodes.items[grid_index].energy = energy;
    grid.nodes.items[grid_index].cache[MB_COORD_X].dirty = true;
    grid.nodes.items[grid_index].cache[MB_COORD_Y].dirty = true;
    grid.nodes.items[grid_index].cache[MB_COORD_Z].dirty = true;
}


double
compute_energy(Vector3d world_pos, BallArray balls) {
    //  e = mass/distance^2 
    
    size_t i;
    double energy = 0;
    for (i = 0; i < balls.count; i++)
    {
        double dist_x = world_pos.x - balls.items[i].pos.x;
        double dist_y = world_pos.y - balls.items[i].pos.y;
        double dist_z = world_pos.z - balls.items[i].pos.z;
        double dist   = (dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);

        if( dist < 0.0001f ) {
            dist = 0.0001f;
        }

        energy += balls.items[i].mass / dist;
        //energy += 1.0f / sqrtf(dist);
    }
    
    return energy;
}

Vertex3d
compute_vertex(Vector3d vertex, Grid grid, BallArray balls) {
    Vertex3d v;
    v.pos = vertex;
    
    // Compute vertex normal
    // n += 2 * mass * vector / distance^4
    Vector3d normal = {0};
    for (size_t i = 0; i < balls.count; i++)
    {
        Vector3d dir = {
            .x = v.pos.x - balls.items[i].pos.x,
            .y = v.pos.y - balls.items[i].pos.y,
            .z = v.pos.z - balls.items[i].pos.z
        };

        double dist = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
        
        normal.x += 2 * balls.items[i].mass * dir.x / (dist * dist);
        normal.y += 2 * balls.items[i].mass * dir.y / (dist * dist);
        normal.z += 2 * balls.items[i].mass * dir.z / (dist * dist);
    }

    // normalize vector
    normal = vector3d_normalize(normal);

    // Calculate the spherical coordinates
    double theta = atan2(normal.z, normal.x);
    double phi   = acos(normal.y);

    // Convert spherical coordinates to texture coordinates
    v.text_coord.x = (theta + PI) / (2.0 * PI);
    v.text_coord.y = phi / PI;
    
    return v;
}

Triangle
compute_triangle(Vector3d vertices[3], Grid grid, BallArray balls) {
    Triangle triangle;

    triangle.vertex1 = compute_vertex(vertices[0], grid, balls);
    triangle.vertex2 = compute_vertex(vertices[1], grid, balls);
    triangle.vertex3 = compute_vertex(vertices[2], grid, balls);

    return triangle;
}

Voxel
compute_voxel(Index3d grid_pos, Grid grid, BallArray balls) {

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

    // gather voxel corner energy values

    Voxel voxel;

    voxel.triangles.count = 0;
    voxel.bits = 0;
    for(size_t c = 0; c < 8; ++c) {
        Index3d pos = {
                        .x = grid_pos.x + marching_cube_vertices[c][0],     
                        .y = grid_pos.y + marching_cube_vertices[c][1],
                        .z = grid_pos.z + marching_cube_vertices[c][2]
                    };
        size_t index = compute_grid_point_index(pos, grid);

        voxel.corners[c].grid_pos   = pos;
        voxel.corners[c].grid_index = index;
        voxel.corners[c].energy     = grid.nodes.items[index].energy;

        voxel.bits |= voxel.corners[c].energy > THRESHOLD ? (1 << c) : 0;
    }

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
// compute edge points to build triangles

    // This table indicates what axis variations have the second index of marching_cube_edges
    size_t coords[12] = {MB_COORD_X, MB_COORD_Z, MB_COORD_X, MB_COORD_Z, 
                        MB_COORD_X, MB_COORD_Z, MB_COORD_X, MB_COORD_Z, 
                        MB_COORD_Y, MB_COORD_Y,  MB_COORD_Y,  MB_COORD_Y};

    Vector3d world_pos = compute_grid_point_coordinates(grid_pos, grid);
    
    Vector3d vertices[3];
    size_t i = 0;
    while(1) {
        size_t vertex_index = i % 3;
        char edge = marching_cube_triangles[voxel.bits][i++];
        if(edge == (char)-1) {
            break;
        }


        char index0 = marching_cube_edges[edge][0];
        char index1 = marching_cube_edges[edge][1];
        Vector3d computedPoint;

#define USE_CACHE
#ifdef USE_CACHE
        size_t coord = coords[edge];
        if(!grid.nodes.items[voxel.corners[index0].grid_index].cache[coord].dirty) {
            computedPoint = grid.nodes.items[voxel.corners[index0].grid_index].cache[coord].value;
        } else {
#endif
        double t = (THRESHOLD - voxel.corners[index0].energy) / (voxel.corners[index1].energy - voxel.corners[index0].energy);
        
        double field_value[MB_COORD_COUNT];
        for(size_t i = 0; i < MB_COORD_COUNT; ++i) {
            double field_intensity_0 = marching_cube_vertices[index0][i];
            double field_intensity_1 = marching_cube_vertices[index1][i];
            field_value[i] = field_intensity_0 + (field_intensity_1 - field_intensity_0) * t;
        }

            computedPoint = (Vector3d)   {
                                    world_pos.x + field_value[MB_COORD_X] * grid.size.x,
                                    world_pos.y + field_value[MB_COORD_Y] * grid.size.y,
                                    world_pos.z + field_value[MB_COORD_Z] * grid.size.z
                                };

#ifdef USE_CACHE
            grid.nodes.items[voxel.corners[index0].grid_index].cache[coord].value = computedPoint;
            grid.nodes.items[voxel.corners[index0].grid_index].cache[coord].dirty = false;
        }
#endif
        vertices[vertex_index] = computedPoint;

        if(i % 3 == 0) {

            Triangle triangle;
            triangle.vertex1 = compute_vertex(vertices[0], grid, balls);
            triangle.vertex2 = compute_vertex(vertices[1], grid, balls);
            triangle.vertex3 = compute_vertex(vertices[2], grid, balls);

            voxel.triangles.items[voxel.triangles.count] = triangle;
            voxel.triangles.count++;
        }
    }
    
    return voxel;
}

// PUBLIC API

void
compute_metaball(BallArray balls, Grid grid) {
    size_t x, y, z;

    for (x = 0; x < grid.count.x; x++) {
        for (y = 0; y < grid.count.y; y++) {
            for (z = 0; z < grid.count.z; z++) {
                Index3d grid_pos = {.x = x, .y = y, .z = z};
                Vector3d world_pos = compute_grid_point_coordinates(grid_pos, grid);
                double energy = compute_energy(world_pos, balls);
                set_energy_point(grid_pos, grid, energy);
            }
        }
    }
}

void
generate_mesh(Camera camera, Grid grid, BallArray balls) {
    size_t x, y, z, i;
    for (x = 0; x < grid.count.x - 1; x++) {
        for (y = 0; y < grid.count.y - 1; y++) {
            for (z = 0; z < grid.count.z - 1; z++) {
                Voxel voxel = compute_voxel((Index3d){.x = x, .y = y, .z = z}, grid, balls);

                for (i = 0; i < voxel.triangles.count; i++) {
                    Vector3 v1 = {voxel.triangles.items[i].vertex1.pos.x, voxel.triangles.items[i].vertex1.pos.y, voxel.triangles.items[i].vertex1.pos.z};
                    Vector3 v2 = {voxel.triangles.items[i].vertex2.pos.x, voxel.triangles.items[i].vertex2.pos.y, voxel.triangles.items[i].vertex2.pos.z};
                    Vector3 v3 = {voxel.triangles.items[i].vertex3.pos.x, voxel.triangles.items[i].vertex3.pos.y, voxel.triangles.items[i].vertex3.pos.z};

                    // Calcular los vectores AB y AC
                    Vector3d AB = {v2.x - v1.x, v2.y - v1.y, v2.z - v1.z};
                    Vector3d AC = {v3.x - v1.x, v3.y - v1.y, v3.z - v1.z};

                    // Calcular la normal del triángulo
                    Vector3d face_normal = vector3d_cross_product(AB, AC);
                    face_normal = vector3d_normalize(face_normal);

                    // vector del centro del triangulo a la camara
                    Vector3d direction = {
                        .x = camera.position.x - (v1.x + v2.x + v3.x) / 3.0,
                        .y = camera.position.y - (v1.y + v2.y + v3.y) / 3.0,
                        .z = camera.position.z - (v1.z + v2.z + v3.z) / 3.0
                    };

                    direction = vector3d_normalize(direction);

                    double dotProduct = vector3d_dot_product(face_normal, direction);

                    DrawTriangle3D(v1, v2, v3, ColorBrightness(RED, dotProduct));
                    DrawLine3D(v1, v2, BLACK);
                    DrawLine3D(v2, v3, BLACK);
                    DrawLine3D(v3, v1, BLACK);

                }
            }
        }
    }
}


void
render(Grid grid, BallArray balls) {
    size_t x, y, z, i;
    
    Vector3 rotationAxis = {.x = 0, .y = 0, .z = 0};

    for (x = 0; x < grid.count.x; x++) {
        for (y = 0; y < grid.count.y; y++) {
            for (z = 0; z < grid.count.z; z++) {
                Vector3 center = {
                                    .x = x * grid.size.x + grid.offset.x, 
                                    .y = y * grid.size.y + grid.offset.y, 
                                    .z = z * grid.size.z + grid.offset.z
                };

                size_t grid_index = compute_grid_point_index((Index3d){.x = x, .y = y, .z = z}, grid);
                if(grid.nodes.items[grid_index].energy <= THRESHOLD) {
                    DrawCircle3D(center, grid.size.x * 0.05f, rotationAxis, 0, BLUE);
                }
                // printf("node[%lu] = [%lu][%lu][%lu] => (%0.2f, %0.2f, %0.2f) -> %0.2f\n", grid_index, x, y, z, center.x, center.y, center.z, grid.nodes.items[grid_index].energy);
            }
        }
    }
    for(i = 0; i < balls.count; ++i) {
        
        Vector3 pos = {
                            .x = balls.items[i].pos.x, 
                            .y = balls.items[i].pos.y, 
                            .z = balls.items[i].pos.z
        };
    }
}

Grid
make_grid(Index3d count, Vector3d size, Vector3d pos) {
    Grid grid;

    Vector3d offset = {
                        .x = pos.x - size.x * ((count.x - 0.5) / 2.0),
                        .y = pos.y - size.y * ((count.y - 0.5) / 2.0),
                        .z = pos.z - size.z * ((count.z - 0.5) / 2.0)
    };

    grid.count  = count;
    grid.size   = size;
    grid.pos    = pos;
    grid.offset = offset;
    grid.nodes  = MB_CREATE_ARRAY(Node, grid.count.x * grid.count.y * grid.count.z);

    return grid;
}

void
free_grid(Grid grid) {

    grid.count  = (Index3d){.x=0, .y=0, .z=0};
    grid.size   = (Vector3d){.x=0, .y=0, .z=0};
    grid.pos    = (Vector3d){.x=0, .y=0, .z=0};
    grid.offset = (Vector3d){.x=0, .y=0, .z=0};
    
    MB_DESTROY_ARRAY(grid.nodes);
}

void
push_ball(BallArray *balls, Ball ball) {
    balls->items[balls->count] = ball;
    balls->count++;
}

void
manage_keys(Grid grid, Camera *camera) {
    static Vector3d movement = {2.5, 0, 1.8};
    double speed = 0.02;
    if(IsKeyDown(KEY_W)) {
        movement.y += speed;
    } else if(IsKeyDown(KEY_S)) {
        movement.y -= speed;
    }
    if(IsKeyDown(KEY_A)) {
        movement.x += speed;
    } else if(IsKeyDown(KEY_D)) {
        movement.x -= speed;
    }
    if(IsKeyDown(KEY_Q)) {
        movement.z += speed;
    } else if(IsKeyDown(KEY_E)) {
        movement.z -= speed;
    }


    camera->position.x = grid.pos.x + sinf(0.5f * movement.x) * ((grid.count.x / 2) * grid.size.x * movement.z);
    camera->position.z = grid.pos.z + sinf(0.7f * movement.x) * ((grid.count.z / 2) * grid.size.z * movement.z);
    camera->position.y = grid.pos.y + sinf(0.5f * movement.y) * ((grid.count.y / 2) * grid.size.y * movement.z);


    // printf("Movement %0.2f, %0.2f, %0.2f \n", movement.x, movement.y, movement.z);
}

void
move_balls(Grid grid, BallArray balls, double t) {
    size_t i = 0;
    {
        balls.items[i].pos.x = grid.pos.x + sinf(0.5f * t) * grid.size.x * 2;
        balls.items[i].pos.y = grid.pos.y + sinf(0.7f * t) * (grid.count.y / 2) * grid.size.y;
        balls.items[i].pos.z = grid.pos.z + sinf(0.2f * t) * (grid.count.z / 4) * grid.size.z;
    }

    i = 2;
    balls.items[i].pos.x =  20 + sinf(2.50f * t) * grid.size.x;
    balls.items[i].pos.y = -20 + sinf(1.75f * t) * grid.size.y;
    balls.items[i].pos.z =  60 + sinf(0.82f * t) * grid.size.z;

}

double
get_current_time() {
#ifdef _WIN32
    LARGE_INTEGER frequency;
    if (!QueryPerformanceFrequency(&frequency)) {
        fprintf(stderr, "ERROR: QueryPerformanceFrequency failed.\n");
        exit(1);
    }

    LARGE_INTEGER now;
    if (!QueryPerformanceCounter(&now)) {
        fprintf(stderr, "ERROR: QueryPerformanceCounter failed.\n");
        exit(1);
    }

    return (double)now.QuadPart / (double)frequency.QuadPart;

#else 
    struct timespec now;
    if (clock_gettime(CLOCK_MONOTONIC, &now) < 0) {
        fprintf(stderr, "ERROR: could not get current monotonic time: %s\n", strerror(errno));
        exit(1);
    }
    return (double) now.tv_sec + now.tv_nsec * 0.000000001f;

#endif
}

int
main(int argc, char *argv[]) {

    Index3d  count = {.x=20, .y=20, .z=20};
    Vector3d size  = {.x=10, .y=10, .z=10};
    Vector3d pos   = {.x=0, .y=0, .z=100};
    Grid    grid  = make_grid(count, size, pos);
    
    BallArray balls = MB_CREATE_ARRAY(Ball, 5);
    push_ball(&balls, (Ball){.pos = {.x=-20, .y=-20, .z=120}, .mass=120});
    push_ball(&balls, (Ball){.pos = {.x=20, .y=20, .z=100}, .mass=300});
    push_ball(&balls, (Ball){.pos = {.x=20, .y=-20, .z=60}, .mass=150});


    const int screenWidth = 800;
    const int screenHeight = 450;

    // Define the camera to look into our 3d world
    Camera camera = { 0 };
    camera.position = (Vector3){ 0.0f, 0.0f, 0.0f }; // Camera position
    camera.target = (Vector3){ 0.0f, 0.0f, 100.0f };      // Camera looking at point
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera projection type

    InitWindow(screenWidth, screenHeight, "Metaballs");
    while (!WindowShouldClose())        // Detect window close button or ESC key
    {

        manage_keys(grid, &camera);
        move_balls(grid, balls, get_current_time());

        compute_metaball(balls, grid);

        BeginDrawing();
            ClearBackground(RAYWHITE);
            BeginMode3D(camera);
            generate_mesh(camera, grid, balls);
            render(grid, balls);
            EndMode3D();
        EndDrawing();
    }
    CloseWindow();

    MB_DESTROY_ARRAY(balls);
    
    free_grid(grid);
    return 0;
}
