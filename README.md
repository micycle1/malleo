# malleo
_As-rigid-as-possible shape deformation for Java_

This library provides a compact, production‑oriented implementation of 2D As‑Rigid‑As‑Possible (ARAP) shape deformation for Java. It implements the method from Igarashi & Igarashi[^1], enabling fast, handle‑based interactive edits of 2D triangle meshes with strong local rigidity preservation.

### Quick example (conceptual)
- Register a mesh (triangulation of your 2D shape).
- Add handles (2D coordinates that act as anchor/manipulation points on the shape)
- Compile/prepare the solver for that handle set.
- On each frame, update the solver with new handle positions and read back deformed shape.

[^1]: Igarashi, T., & Igarashi, Y. (2009). Implementing As‑Rigid‑As‑Possible Shape Manipulation and Surface Flattening. Journal of Graphics, GPU, and Game Tools, 14(1), 17–30. https://doi.org/10.1080/2151237X.2009.10129273
