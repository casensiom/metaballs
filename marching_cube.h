#pragma once
#ifndef _MARCHINGCUBE_H_
#define _MARCHINGCUBE_H_

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
static float marching_cube_vertices[8][3] = 
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




#endif