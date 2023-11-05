# Metaballs

![Thumbnail](./thumbnail.jpg)

Header only 3d metaballs generator in ANSI C

## About

[![Build Status](https://github.com/casensiom/metaballs/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/casensiom/metaballs/actions)

This library computes 3d metaballs in pure ANSI C.

![Thumbnail](./thumbnail.gif)

There are optimizations yet to be done.

The example uses [raylib](https://www.raylib.com) as a render engine.

## How to use

This library follows [single-file C/C++ public domain libraries](https://github.com/nothings/single_file_libs) rules.

Include the header file and the definition to include the implementation code in just one of the includes.

```c
#define METABALLS_IMPLEMENTATION
#include "metaballs.h"
``` 

Refer [main.c](https://github.com/casensiom/metaballs/main.c) to see the library usage.

## Features
- **Iteration**: Added a method to visit only valid neighbors.
- **Optimization**: A rudimentary cache system is included.
- **Color**: Added color per ball.
- **Texture mapping**: Generate and upload mesh to GPU, texture mapping is available.
- **Control**: Added mouse and keyboard control.

## TODO
 - Code clean: Add prefix to library methods and structures. (Avoid name collision)
 - Known issue: The spheric texture mapping doesn't close correctly.
 - Performance: Avoid repeat vertex on mesh generation.
