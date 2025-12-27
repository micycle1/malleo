[![](https://jitpack.io/v/micycle1/malleo.svg)](https://jitpack.io/#micycle1/malleo)

# malleo
_As-rigid-as-possible shape deformation for Java_

This library provides a compact, production‑oriented implementation of 2D As‑Rigid‑As‑Possible (ARAP) shape deformation for Java. It implements the method from Igarashi & Igarashi[^1], enabling fast, handle‑based interactive edits of 2D triangle meshes with strong local rigidity preservation.

### Quick example (conceptual)
- Register a mesh (triangulation of your 2D shape).
- Add handles (2D coordinates that act as anchor/manipulation points on the shape)
- Compile/prepare the solver for that handle set.
- On each frame, update the solver with new handle positions and read back deformed shape.

### Quick Start

Malleo is built to interop with JTS Geometries 
```Java
// 1) Rest shape (square)
Geometry rest = new WKTReader().read(
  "POLYGON ((200 200, 600 200, 600 600, 200 600, 200 200))"
);

// 2) Triangulate
Geometry triangles = ConstrainedDelaunayTriangulator.triangulate(rest);

// 3) Create deformer
Malleo malleo = new Malleo(triangles);

// 4) Define handles (rest-space). Order matters.
List<Coordinate> handles = new ArrayList<>(List.of(
  new Coordinate(200, 200),
  new Coordinate(600, 600)  // (will move)
));

// 5) Compile once (in a real sketch, keep this + malleo as fields and reuse)
Malleo.CompiledHandles compiled = malleo.prepareHandles(handles);

// 6) Move ONLY the second handle and solve
handles.get(1).setX(700);
handles.get(1).setY(200);

Polygon deformed = malleo.solve(compiled, handles);
```

[^1]: Igarashi, T., & Igarashi, Y. (2009). Implementing As‑Rigid‑As‑Possible Shape Manipulation and Surface Flattening. Journal of Graphics, GPU, and Game Tools, 14(1), 17–30. https://doi.org/10.1080/2151237X.2009.10129273
