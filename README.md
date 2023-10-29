# Metaballs
![Thumbnail](./thumbnail.png)

## Description

This proyect builds a metaball in pure C, uses raylib to render it.

There are optimizations to be done.

## TODO

 - Iteration, currently we iterate over each grid voxel, improve it.
    > We need to find one of the voxels at the energy field limit, then compute the voxel and iterate over their neighbors. Repeat for each ball.

 - ~~Compute the edge values only once~~. A rudimentary cache system is included.
 - Clean code, separate voxel calculation and render code.
 - Control, add mouse and keyboard control.
