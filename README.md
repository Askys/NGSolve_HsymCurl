An implementation of novel H(symCurl) conforming finite elements in NGSolve.

To use, add the files from one of the sub-folders (Y0,S0 or Y1) to your main folder.

An example (also in available in conv_benchmark.ipynb):

```python
# Compile the module
from ngsolve.fem import CompilePythonModule
from pathlib import Path

txt = Path('mymodule.cpp').read_text()
m = CompilePythonModule(Path('mymodule.cpp'), init_function_name='mymodule')

# Build a mesh
refinement = 1
l = refinement * 4
mesh = MakeStructured3DMesh(False, nx=3 * int(l), ny=int(1 * l), nz=int(1 * l),
                            mapping=lambda x, y, z: (-3 + 6 * x, -1 + 2 * y, -1 + 2 * z))
Draw(mesh)

# Define the finite element space
fesHsC = m.MyNedMatFESpace(mesh, order=1)  # for S1, otherwise m.HsymCurlFESpace(mesh, order=0) for Y0 or S0 
fes = fesHsC

print(fes.ndof)

# Define the force
tol = 1e-10

# Jumping constant identity
O = CF((0, 0, 0,
        0, 0, 0,
        0, 0, 0), dims=(3, 3))

I = CF((1, 0, 0,
        0, 1, 0,
        0, 0, 1), dims=(3, 3))

R = IfPos(Norm(x) - 1 - tol, O, I)

# Jumping idenity
# O = CF( (0, 0, 0,
#          0, 0, 0,
#          0, 0, 0), dims=(3,3) )

# I = CF( (1, 0, 0,
#          0, 1, 0,
#          0, 0, 1), dims=(3,3) )

# V = sin(x + 2. * y - 3. * z) * I

# R = IfPos( Norm(x) - 1 - tol, V, 2*V  )

# Jumping normal
# O = CF( (0, 0, 0,
#          0, 0, 0,
#          0, 0, 0), dims=(3,3) )

# I = CF( (1, 0, 0,
#          1, 0, 0,
#          1, 0, 0), dims=(3,3) )

# R = IfPos( Norm(x) - 1 - tol, sin(x)*I ,  cos(x)*I  )

# Contininous field
# O = CF( (0, 0, 0,
#          0, 0, 0,
#          0, 0, 0), dims=(3,3) )

# I = CF( (1, 0, 0,
#          0, 0, 0,
#          0, 0, 0), dims=(3,3) )

# R = sinh(x)/10*I

clipping = {"pnt": (0, 0, 0), "vec": (0, 1, 0)}
Draw(R, mesh, clipping=clipping)

# Construct a bilinear form
P, dP = fes.TnT()
a = BilinearForm(fes, symmetric=True, symmetric_storage=True, condense=False)
a += (InnerProduct(dP, P)
      + InnerProduct(curl(dP), curl(P))  # Actually symCurl
      ) * dx

f = LinearForm(fes)
f += (InnerProduct(dP, R)) * dx

sol = GridFunction(fes)

# Solve
sol.vec[:] = 0

r = sol.vec.CreateVector()
w = sol.vec.CreateVector()

with TaskManager():
    f.Assemble()
    a.Assemble()
    inv = a.mat.Inverse(fes.FreeDofs(a.condense), inverse="sparsecholesky")

    a.Apply(sol.vec, r)

    r.data -= f.vec
    if a.condense:
        r.data += a.harmonic_extension_trans * r
    w.data = inv * r
    if a.condense:
        w.data += a.harmonic_extension * w
        w.data += a.inner_solve * r
    sol.vec.data -= w

# Draw the solution 
gfR0 = CF((sol[0, 0], sol[0, 1], sol[0, 2]))
gfR1 = CF((sol[1, 0], sol[1, 1], sol[1, 2]))
gfR2 = CF((sol[2, 0], sol[2, 1], sol[2, 2]))

clipping = {"pnt": (0, 0, 0), "vec": (0, 1, 0)}
Draw(gfR0, mesh, "P", clipping=clipping, draw_surf=False, vectors={"grid_size": 100})
Draw(gfR1, mesh, "P", clipping=clipping, draw_surf=False, vectors={"grid_size": 100})
Draw(gfR2, mesh, "P", clipping=clipping, draw_surf=False, vectors={"grid_size": 100})
```