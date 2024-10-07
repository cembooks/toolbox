<h1>Test program for observing shape functions</h1>

<h2> Introduction </h2>

This program was used to compose the swap tables in the pull "Nedelec orientation
fix nr.1". It can also be used for verification. 

<h2> Operation </h2>

For this program to work properly, one needs to comment out the content of a 
function in grid_tools_topology.cc:

@code
  template <int dim>
  void
  consistently_order_cells(std::vector<CellData<dim>> &cells)
  {
/*
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
*/  
  }
@endcode

![][fig-shape-finctions]

[fig-shape-finctions]: doc/shape-functions.png

