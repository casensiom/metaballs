#Metaballs
---
![Thumbnail](./thumbnail.png)

# Description

This proyect builds a metaball in pure C, uses raylib to render it.

There are optimizations to be done.

## TODO

 - Compute the edge values only once. We are now computing all the edges of each voxel, allowing contiguous voxels to reuse previous calculations.
 - Clean code, separate voxel calculation and render code.