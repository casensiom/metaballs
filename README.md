# Metaballs

Header only metaballs generator in ANSI C

![Thumbnail](./thumbnail.gif)

## About

[![Build Status](https://github.com/casensiom/metaballs/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/casensiom/metaballs/actions)


This proyect computes metaballs in pure ANSI C.

There are optimizations yet to be done.

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
