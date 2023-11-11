#ifndef _BALL_H_
#define _BALL_H_

#include "../metaballs.h"

void ball_push(BallArray *balls, Ball ball);
Ball ball_create_random(Grid grid);

#endif //_BALL_H_