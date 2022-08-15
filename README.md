
# vector graph, delaunay triangulation, voronoi
todo

## Features
* Vector graph
* Delaunay triangulation with a parallel algorithm
* Build voronoi graph
* todo ...

## Usage
`todo.cpp`

```CPP
#include "vector_graph.h"
#include "vector_graph_delaunay.h"
#include "voronoi_graph.h"

#include <array>

// custom vector 2D type with operator[](int)
struct Vec2d : public std::array<double, 2>
{
    Vec2d(double X, double Y) : std::array<double, 2>({X, Y}) {}
    Vec2d(double In) : std::array<double, 2>({In, In}) {}
}

// example input verts
std::vector<Vec2d> Verts = GetVerts2d();

// build delaunay triangulated graph
auto TriangulatedGraph = Triangulate::triangulate<Vec2d>(Verts);

// build voronoi graph
auto VoronoiGraph = Voronoi::voronoi<Vec2d>(TriangulatedGraph);

//draw graph example:
TriangulatedGraph.for_each_edge(
    [&](const Vec2d& V1, const Vec2d& V2)
    {
        // some DrawLine API for draw graph
        DrawLine(V1, V2 /*, color etc*/);
    });
//draw graph example:
VoronoiGraph.for_each_edge(
    [&](const Vec2d& V1, const Vec2d& V2)
    {
        // some DrawLine API for draw graph
        DrawLine(V1, V2 /*, color etc*/);
    });


/********/


//support custom vector 3D type with operator[](int)
struct Vec3d : public std::array<double, 3>
{
    Vec3d(double X, double Y, double Z) : std::array<double, 3>({X, Y, Z}) {}
    Vec3d(double In) : std::array<double, 3>({In, In, In}) {}
}
std::vector<Vec3d> Verts3d = GetVerts3d();
// flat z triangulation 
// based on Vec3d : X -[0] , y - [1], z - [2] ignored
auto TriangulatedGraph3d = Triangulate::triangulate<Vec3d>(Verts3d);
auto VoronoiGraph3d = Voronoi::voronoi<Vec3d>(TriangulatedGraph3d);

```

## Benchmarks
```
todo
```