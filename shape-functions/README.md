<h1>Test program for observing shape functions</h1>

<h2> Introduction </h2>

This program was used to compose the swap tables in the pull "Nedelec orientation
fix nr.1". It can also be used for verification.

<h2> Operation </h2>

For this program to work properly, one needs to comment out the content of a
function in 
[grid_tools_topology.cc](https://github.com/dealii/dealii/blob/master/source/grid/grid_tools_topology.cc#L1613):

    template <int dim>
    void
    consistently_order_cells(std::vector<CellData<dim>> &cells)
    {/*
       Assert(cells.size() != 0,
         ExcMessage(
             "List of elements to orient must have at least one cell"));

       // there is nothing for us to do in 1d
       if (dim == 1)
         return;

       // check if grids are already consistent. if so, do
       // nothing. if not, then do the reordering
       if (!is_consistent(cells))
         try
         {
           reorient(cells);
         }
         catch (const ExcMeshNotOrientable &)
         {
          // the mesh is not orientable. this is acceptable if we are in 3d,
          // as class Triangulation knows how to handle this, but it is
          // not in 2d; in that case, re-throw the exception
          if (dim < 3)
            throw;
          }
     */}

This function being enabled reorders the mesh what makes test of certain face orientations impossible.

The shape-functions program allows to plot shape functions on three meshes: two separate one-cell meshes
and one two-cells mesh. The two-cells mesh is derived by merging together the two separate one-cell meshes.
The figure below illustrates such plots.

![][fig-shape-finctions]

The plot on the left illustrates an overlay of the first shape functions on the shared face of the two
separate one-cell meshes. From this plot one can deduce that a swap of the shape functions is needed.
The plot on the right illustrates the corresponding shape function on the two-cells mesh. It suggests
that the swap and the sign change encoded in the swap table are correct as the shape functions in both
cells are oriented consistently.

All controls of the shape-functions program are located in the last two lines of the 
[CMakeLists.txt](https://github.com/cembooks/toolbox/blob/main/shape-functions/CMakeLists.txt):

    target_compile_options(${TARGET} PRIVATE -DDIMENSION__=2 -DFEDEGREE__=2
         -DFACEORIENTATION__=0)

The macro definition DIMENSION__ can take two values: 2 and 3.  It corresponds to the parameter dim
in deal.II.

The macro definition FEDEGREE__ must be non-negative integer. This is the degree of the Nedelec
finite element being tested, i.e., FE_NEDELEC<DIMENSION__> fe(FEDEGREE__).

The macro definition FACEORIENTATION__ in two dimensions encodes the orientation of the shared edge:
0 - opposite orientation, 1 - normal orientation. In three dimensions the macro definition
FACEORIENTATION__ can take the following values: 0, 1, 2, 3, 4, 5, 6, or 7. This numbers are
the combined orientation and are computed as

    orientation_no = face_flip*4 + face_rotation*2 + face_orientation*1;

See
[tria_orientation.h](https://github.com/dealii/dealii/blob/master/include/deal.II/grid/tria_orientation.h).
The three boolean parameters in the last equation have the following meaning:

>face_orientation - true if face is in standard orientation and false otherwise;  
>face_rotation - rotation by 90 deg counterclockwise if true;  
>face_flip - rotation by 180 deg counterclockwise if true.  

See the documentation of
[GeometryInfo<dim>](https://dealii.org/current/doxygen/deal.II/structGeometryInfo.html).



[fig-shape-finctions]: doc/shape-functions.png

