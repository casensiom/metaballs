#include "ball.h"

void ball_push(BallArray *balls, Ball ball)
{
    balls->items[balls->count] = ball;
    balls->count++;
}

Ball ball_create_random(Grid grid)
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

    return (Ball){.pos = p, .mass = ball_mass, .color = color};
}
