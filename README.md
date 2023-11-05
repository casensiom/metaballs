# Metaballs
![Thumbnail](./thumbnail.gif)

## Description

[![Build Status](https://github.com/casensiom/metaballs/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/casensiom/metaballs/actions)


This proyect computes metaballs in pure C, using raylib to render it.

There are optimizations to be done.

## TODO

 - ~~Iteration~~. Added a method to visit only valid neighbors.
 - ~~Compute the edge values only once~~. A rudimentary cache system is included.
 - ~~Add color~~. Added color per ball.
 - ~~Add texture mapping~~. Generate and upload mesh to GPU, texture mapping is available.
 - ~~Control~~. Added mouse and keyboard control.
 - Clean code, separate voxel calculation and render code.
