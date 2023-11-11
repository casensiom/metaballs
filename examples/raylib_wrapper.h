
#include "../metaballs.h"

#include <raylib.h>
#include <rlgl.h>


Mesh raylib_upload_mesh(Grid grid);
void raylib_render_mesh(Grid grid, Mesh *mesh, Material material, Material wired, Matrix transform);
void raylib_unload_mesh(Mesh *mesh);
