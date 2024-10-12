/******************************************************************************
 * Copyright (C) Siarhei Uzunbajakau, 2023.
 *
 * This program is free software. You can use, modify, and redistribute it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 or (at your option) any later version.
 * This program is distributed without any warranty.
 *
 * Refer to COPYING.LESSER for more details.
 ******************************************************************************/

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_orientation.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace dealii;

template <int dim>
class ShapeFunctions
{
public:
  ShapeFunctions();

  void
  run();

private:
  const unsigned int combined_face_orientation = FACEORIENTATION__;

  void
  make_mesh();
  
  void
  print_shared_face(TriaActiveIterator<CellAccessor<dim, dim>> cell);
  
  void
  setup_system();
  
  void
  save(std::string fname) const;

  Triangulation<dim> triangulation;
  FE_Nedelec<dim>    fe;
  DoFHandler<dim>    dof_handler;
  SparsityPattern    sparsity_pattern;
  Vector<double>     solution;

  const std::string fname_vtk =
    (DIMENSION__ == 2) ? "Data/2D_shape_function" : "Data/3D_shape_function";
};

template <int dim>
ShapeFunctions<dim>::ShapeFunctions()
  : fe(FEDEGREE__)
{
  if (dim == 2)
    if (combined_face_orientation > 3)
      {
        std::cout << "Error. In the two-dimensional space the macro definition \
FACEORIENTATION__ must be in the range 0...3.\n";
        return;
      }

  if (dim == 3)
    if (combined_face_orientation > 7)
      {
        std::cout
          << "Error. In the three-dimensional space the macro definition \
FACEORIENTATION__ must be in the range 0...7.\n";
        return;
      }

  run();
}

template <>
void
ShapeFunctions<2>::print_shared_face(
  TriaActiveIterator<CellAccessor<2, 2>> cell)
{
  std::cout << "Cell id = " << cell->id() << "-----------------------------"
            << std::endl;

  for (unsigned int f = 0; f < 4; f++)
    if (cell->neighbor_index(f) != -1)
      {
        std::cout << "The local index of the shared face: " << f << std::endl;
        std::cout << "Shared face - orientation: " << cell->line_orientation(f)
                  << std::endl;
      }
}

template <>
void
ShapeFunctions<3>::print_shared_face(
  TriaActiveIterator<CellAccessor<3, 3>> cell)
{
  unsigned int face = -1;

  for (unsigned int f = 0; f < 6; f++)
    if (cell->neighbor_index(f) != -1)
      face = f;

  std::cout << "Cell id = " << cell->id() << ": " << std::endl;

  std::cout << "The local index of the shared face: " << face << std::endl;

  std::cout << "Shared face - flip: " << cell->face_flip(face) << std::endl;
  std::cout << "Shared face - rotation: " << cell->face_rotation(face)
            << std::endl;
  std::cout << "Shared face - orientation: " << cell->face_orientation(face)
            << std::endl;
  std::cout << "Shared face - combined orientation: "
            << static_cast<unsigned int>(internal::combined_face_orientation(
                 cell->face_orientation(face),
                 cell->face_rotation(face),
                 cell->face_flip(face)))
            << std::endl;
}

template <>
void
ShapeFunctions<2>::make_mesh()
{
  GridGenerator::non_standard_orientation_mesh(triangulation,
                                               combined_face_orientation);

  for (auto cell : triangulation.active_cell_iterators())
    print_shared_face(cell);
}


template <>
void
ShapeFunctions<3>::make_mesh()
{
  bool face_orientation;
  bool face_rotation;
  bool face_flip;

  std::tie(face_orientation, face_rotation, face_flip) =
    dealii::internal::split_face_orientation(combined_face_orientation);

  std::cout << "Orientation ask: " << face_flip << face_rotation
            << face_orientation << std::endl;

  GridGenerator::non_standard_orientation_mesh(
    triangulation, face_orientation, face_flip, face_rotation, false);

  for (auto cell : triangulation.active_cell_iterators())
    print_shared_face(cell);
}

template <int dim>
void
ShapeFunctions<dim>::setup_system()
{
  dof_handler.reinit(triangulation);
  dof_handler.distribute_dofs(fe);

  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  AffineConstraints<double>(),
                                  false);

  sparsity_pattern.copy_from(dsp);
  solution.reinit(dof_handler.n_dofs());
}

template <int dim>
void
ShapeFunctions<dim>::save(std::string fname) const
{
  std::vector<std::string> solution_names(dim, "ShapeFunction");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim,
                   DataComponentInterpretation::component_is_part_of_vector);

  DataOut<dim> data_out;

  data_out.add_data_vector(dof_handler,
                           solution,
                           solution_names,
                           interpretation);

  data_out.build_patches(16);

  std::ofstream out(fname);

  data_out.write_vtk(out);
  out.close();
}

template <int dim>
void
ShapeFunctions<dim>::run()
{
  make_mesh();
  setup_system();

  std::cout << std::endl << fe.get_name() << std::endl;

  std::cout << "\nWriting to ./Data/\n...\n";

  unsigned int N = solution.size();
  for (unsigned int i = 0; i < N; i++)
    {
      solution.reinit(N);
      solution(i) = 1.0;
      save(fname_vtk + std::to_string(i) + ".vtk");
    }
}

int
main()
{
  ShapeFunctions<DIMENSION__> shape_functions;

  return 0;
}
