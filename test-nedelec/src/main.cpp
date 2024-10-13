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

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_orientation.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace dealii;

template <int dim>
class Weight : public Function<dim>
{
public:
  virtual double
  value(const Point<dim>  &r,
        const unsigned int component = 0) const override final;
};

class MainOutputTable : public ConvergenceTable
{
public:
  MainOutputTable() = delete;
  MainOutputTable(int dimensions)
    : ConvergenceTable()
    , dimensions(dimensions)
  {}

  void
  save(std::string fname);

  void
  set_new_order(std::vector<std::string> new_order_in)
  {
    new_order = new_order_in;
  }

  void
  append_new_order(std::string new_clmn)
  {
    new_order.push_back(new_clmn);
  }

  virtual void
  format();

private:
  int                      dimensions;
  std::vector<std::string> new_order =
    {"p", "r", "ncells", "ndofs", "L2"};
};

template <int dim>
class MagneticVectorPotential : public Function<dim>
{
public:
  MagneticVectorPotential()
    : Function<dim>(dim){};

  virtual void
  vector_value_list(const std::vector<Point<dim>> &r,
                    std::vector<Vector<double>>   &values) const override final;

private:
  const double pi =
    3.141592653589793238462643383279502884197169399375105820974944592307816406286;
  const double k = 0.5 * pi;
};

template <int dim>
class TestNedelec
{
public:
  TestNedelec() = delete;
  TestNedelec(unsigned int degree, unsigned int nr_ref);

  void
  run();

  unsigned int
  get_n_cells() const
  {
    return static_cast<unsigned int>(triangulation.n_active_cells());
  }

  unsigned int
  get_n_dofs() const
  {
    return static_cast<unsigned int>(dof_handler.n_dofs());
  }

  double
  get_L2_norm() const
  {
    return L2_norm;
  }

private:
  double        L2_norm;
  Vector<float> L2_per_cell;

  const unsigned int combined_face_orientation = FACEORIENTATION__;
  const unsigned int number_of_mesh_refinements;

  Triangulation<dim> triangulation_right;
  Triangulation<dim> triangulation_left;
  Triangulation<dim> triangulation;

  FE_Nedelec<dim>              fe;
  MagneticVectorPotential<dim> magnetic_vector_potential;

  DoFHandler<dim> dof_handler;

  SparsityPattern           sparsity_pattern;
  SparseMatrix<double>      system_matrix;
  AffineConstraints<double> constraints;

  const FEValuesExtractors::Vector VE;

  Vector<double> solution;
  Vector<double> system_rhs;

  const std::string fname_vtk = "Data/projection";

  void
  make_mesh();

  void
  setup_system();

  void
  assemble_system();

  void
  solve();

  void
  compute_error_norms();

  void
  save() const;
};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
template <int dim>
double
Weight<dim>::value(const Point<dim> &r, const unsigned int component) const
{
  return 1.0;
}
#pragma GCC diagnostic pop

void
MainOutputTable::format()
{
  set_precision("L2", 2);

  set_scientific("L2", true);

  evaluate_convergence_rates("L2",
                             "ncells",
                             ConvergenceTable::reduction_rate_log2,
                             dimensions);

  set_tex_caption("p", "p");
  set_tex_caption("r", "r");
  set_tex_caption("ncells", "nr. cells");
  set_tex_caption("ndofs", "nr. dofs");
  set_tex_caption("L2", "L2 norm");

  set_column_order(new_order);
}

void
MainOutputTable::save(std::string fname)
{
  format();

  std::cout << "------------------------------\n";
  write_text(std::cout);
  std::cout << "\n\n";

  // Save the table in text ...
  {
    std::ofstream ofs(fname + ".txt");
    write_text(ofs);
  }
  // and tex formats.
  {
    std::ofstream ofs(fname + ".tex");
    write_tex(ofs);
  }
}

template <>
void
MagneticVectorPotential<3>::vector_value_list(
  const std::vector<Point<3>> &r,
  std::vector<Vector<double>> &values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  auto v = values.begin();
  for (auto p : r)
    {
      (*v)(0) = -sin(k * p[1]) / k;
      (*v)(1) = sin(k * p[0]) / k;
      (*v)(2) = 0.0;
      v++;
    }
}

template <>
void
MagneticVectorPotential<2>::vector_value_list(
  const std::vector<Point<2>> &r,
  std::vector<Vector<double>> &values) const
{
  Assert(values.size() == r.size(),
         ExcDimensionMismatch(values.size(), r.size()));

  auto v = values.begin();
  for (auto p : r)
    {
      (*v)(0) = -sin(k * p[1]) / k;
      (*v)(1) = sin(k * p[0]) / k;

      v++;
    }
}

template <int dim>
TestNedelec<dim>::TestNedelec(unsigned int degree, unsigned int nr_ref)
  : number_of_mesh_refinements(nr_ref)
  , fe(degree)
  , dof_handler(triangulation)
  , VE(0)
{}

template <>
void
TestNedelec<2>::make_mesh()
{
  GridGenerator::non_standard_orientation_mesh(triangulation,
                                               combined_face_orientation);

  triangulation.refine_global(number_of_mesh_refinements);
}

template <>
void
TestNedelec<3>::make_mesh()
{
  bool face_orientation;
  bool face_rotation;
  bool face_flip;

  std::tie(face_orientation, face_rotation, face_flip) =
    dealii::internal::split_face_orientation(combined_face_orientation);

  GridGenerator::non_standard_orientation_mesh(
    triangulation, face_orientation, face_flip, face_rotation, false);

  triangulation.refine_global(number_of_mesh_refinements);
}

template <int dim>
void
TestNedelec<dim>::setup_system()
{
  dof_handler.reinit(triangulation);
  dof_handler.distribute_dofs(fe);

  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}

template <int dim>
void
TestNedelec<dim>::assemble_system()
{
  QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  std::vector<Vector<double>> exact_solution_values(
    fe_values.get_quadrature().size(), Vector<double>(dim));

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);

      magnetic_vector_potential.vector_value_list(
        fe_values.get_quadrature_points(), exact_solution_values);

      cell_matrix = 0;
      cell_rhs    = 0;

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          for (const unsigned int i : fe_values.dof_indices())
            {
              for (const unsigned int j : fe_values.dof_indices())
                cell_matrix(i, j) += fe_values[VE].value(i, q_index) *
                                     fe_values[VE].value(j, q_index) *
                                     fe_values.JxW(q_index);

              for (unsigned int k = 0; k < dim; k++)
                cell_rhs(i) += fe_values[VE].value(i, q_index)[k] *
                               exact_solution_values.at(q_index)[k] *
                               fe_values.JxW(q_index);
            }
        }
      cell->get_dof_indices(local_dof_indices);

      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
}

template <int dim>
void
TestNedelec<dim>::solve()
{
  ReductionControl control(1000000, 0.0, 1e-12, false, false);

  GrowingVectorMemory<Vector<double>> memory;
  SolverCG<Vector<double>>            cg(control, memory);

  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve(system_matrix, solution, system_rhs, preconditioner);

  constraints.distribute(solution);
}

template <int dim>
void
TestNedelec<dim>::save() const
{
  std::vector<std::string> solution_names(dim, "MagneticVectorPotential");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim,
                   DataComponentInterpretation::component_is_part_of_vector);

  DataOut<dim> data_out;

  data_out.add_data_vector(dof_handler,
                           solution,
                           solution_names,
                           interpretation);

  data_out.build_patches(fe.degree + 2);

  std::ofstream out(fname_vtk + std::to_string(DIMENSION__) + "D_p" +
                    std::to_string(fe.degree - 1) + "_r " +
                    std::to_string(number_of_mesh_refinements) + ".vtk");

  data_out.write_vtk(out);
  out.close();
}

template <int dim>
void
TestNedelec<dim>::compute_error_norms()
{
  Weight<dim>                  weight;
  const Function<dim, double> *mask = &weight;

  QGauss<dim> quadrature(fe.degree + 3);

  VectorTools::integrate_difference(MappingQ<dim>(1),
                                    dof_handler,
                                    solution,
                                    magnetic_vector_potential,
                                    L2_per_cell,
                                    quadrature,
                                    VectorTools::L2_norm,
                                    mask);

  L2_norm = VectorTools::compute_global_error(triangulation,
                                              L2_per_cell,
                                              VectorTools::L2_norm);
}

template <int dim>
void
TestNedelec<dim>::run()
{
  make_mesh();
  setup_system();
  assemble_system();
  solve();
  compute_error_norms();
  save();
}

int
main()
{
  std::cout << "Dimensions: " << DIMENSION__ << std::endl
            << "Face orientation: " << FACEORIENTATION__ << std::endl
            << "FE degree: ";

  std::vector<MainOutputTable> tables(5, MainOutputTable(DIMENSION__));

  for (unsigned int p = 0; p < 5; p++)
    {
      std::cout << p << " ";
      std::flush(std::cout);

      unsigned int r0;
#if DIMENSION__ == 2
      switch (p)
        {
          case 0:
            r0 = 4;
            break;
          case 1:
            r0 = 3;
            break;
          case 2:
            r0 = 3;
            break;
          case 3:
            r0 = 2;
            break;
          case 4:
            r0 = 1;
            break;
        }
#endif

#if DIMENSION__ == 3
      switch (p)
        {
          case 0:
            r0 = 2;
            break;
          case 1:
            r0 = 2;
            break;
          case 2:
            r0 = 1;
            break;
          case 3:
            r0 = 1;
            break;
          case 4:
            r0 = 0;
            break;
        }
#endif

      tables.at(p).clear();
      for (unsigned int r = 0; r < 3; r++)
        {
          tables.at(p).add_value("r", r);
          tables.at(p).add_value("p", p);

          TestNedelec<DIMENSION__> test(p, r0 + r);
          test.run();

          tables.at(p).add_value("ndofs", test.get_n_dofs());
          tables.at(p).add_value("ncells", test.get_n_cells());
          tables.at(p).add_value("L2", test.get_L2_norm());
        }
    }

  std::cout << std::endl;

  for (unsigned int p = 0; p < 5; p++)
    tables.at(p).save("Data/main_table_p" + std::to_string(p));

  return 0;
}
