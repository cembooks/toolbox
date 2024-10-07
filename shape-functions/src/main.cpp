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

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_orientation.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace dealii;

template <int dim>
class ShapeFunctions
{
	public:
		ShapeFunctions();

		void run();

	private:
	
		const unsigned int combined_face_orientation = FACEORIENTATION__;

		void make_mesh();
		void setup_system_left();
		void setup_system_right();
		void setup_system();

		void save_left(std::string fname) const;
		void save_right(std::string fname) const;
		void save(std::string fname) const;

		Triangulation<dim> triangulation_right;
		Triangulation<dim> triangulation_left;
		Triangulation<dim> triangulation;

		FE_Nedelec<dim> fe;

		DoFHandler<dim> dof_handler_left;
		DoFHandler<dim> dof_handler_right;
		DoFHandler<dim> dof_handler;

		SparsityPattern sparsity_pattern_left;
		SparsityPattern sparsity_pattern_right;
		SparsityPattern sparsity_pattern;

		Vector<double> solution_left;
		Vector<double> solution_right;
		Vector<double> solution;

		const std::string fname_vtk = ( DIMENSION__ == 2 ) ?
			"Data/2D_shape_function" : "Data/3D_shape_function";
};

template<int dim>
ShapeFunctions<dim>::ShapeFunctions()
	: fe(FEDEGREE__), 
	dof_handler_left(),
	dof_handler_right()
{
	run();
}

template<>
void ShapeFunctions<2>::make_mesh()
{
	std::cout << "\n\n";
	std::cout << "                                  \n";
 	std::cout << "    2        1        2(4)        \n";   
 	std::cout << "     *-------*-------*            \n";
 	std::cout << "     |       |       |            \n";
 	std::cout << "     |       |       |            \n";
 	std::cout << "     |       |       |            \n";
 	std::cout << "     *-------*-------*            \n";
	std::cout << "    3        0        3(5)        \n";
 	std::cout << "                                  \n";
	std::cout << "                                  \n";
	std::cout << "\n\n";

	std::vector<Point<2>> vertices_left =
	{{ 0.0, 0.0 },  // 0
	 { 0.0, 1.0 },  // 1
	 {-1.0, 1.0 },  // 2
	 {-1.0, 0.0 }}; // 3

	std::vector<Point<2>> vertices_right =
	{{ 0.0, 0.0 },  // 0
	 { 0.0, 1.0 },  // 1
	 { 1.0, 1.0 },  // 2
	 { 1.0, 0.0 }}; // 3

	CellData<2> cell_left, cell_right;

	if (combined_face_orientation % 2 == 0)
	{
		cell_left.vertices = {3, 0, 2, 1};
		cell_right.vertices = {2, 1, 3, 0};
	}
	else
	{
		cell_left.vertices = {3, 0, 2, 1};
		cell_right.vertices = {0, 3, 1, 2};
 }

	std::cout << "Two one-cell meshes --------------------------------\n";

	std::cout << "Left-mesh cell  | ";
	for (unsigned j : {0, 1, 2, 3})
		std::cout << std::setw(3) << cell_left.vertices.at(j);
	std::cout << std::endl;

	std::cout << "Right-mesh cell | ";
	for (unsigned j : {0, 1, 2, 3})
		std::cout << std::setw(3) << cell_right.vertices.at(j);
	std::cout << std::endl;
	std::cout << std::endl;

	std::vector<CellData<2>> cell_data_left = {cell_left};
	std::vector<CellData<2>> cell_data_right = {cell_right};

	SubCellData subcell_data;  

	triangulation_left.create_triangulation(vertices_left, cell_data_left, subcell_data);
	triangulation_right.create_triangulation(vertices_right, cell_data_right, subcell_data);
	GridGenerator::merge_triangulations(triangulation_left, triangulation_right, triangulation);

	// Now we have three meshes:
	// triangulation_left - single cell in the second quadrant; the origin is a vertex. 
	// triangulation_right - single cell in the first quadrant; the origin is a vertex. 
	// triangulation - exactly the same two cells combined in one mesh. 

	// Save the combined mesh.
	GridOut gridout;
	GridOutFlags::Msh msh_flags(true, true);
	gridout.set_flags(msh_flags);

	std::ofstream ofs("output.msh");
	gridout.write_msh(triangulation, ofs);

  // Let's print diagnostic information on the two-cells mesh, triangulation.
	
	std::cout << "Two-cells mesh -------------------------------------\n";
	
	std::tuple< std::vector<Point<2>>, std::vector<CellData<2>>, SubCellData> mesh_description;
	mesh_description = GridTools::get_coarse_mesh_description(triangulation);

	std::cout << "Vertices:\n";
	for (unsigned int i = 0; i < std::get<0>(mesh_description).size(); i++)
	{
		std::cout << std::setw(4) << i << " | ";

		for (unsigned j : {0, 1})
			std::cout << std::setw(4) << std::right << std::get<0>(mesh_description).at(i)[j];

		std::cout << std::endl;
	}

	std::cout << "Cells:\n";
	for (unsigned int i : {0, 1})
	{
		std::cout << "Cell id = " << i << " | ";
	
		for (unsigned int j = 0; j < 4; j++)
			std::cout << std::setw(3) << std::get<1>(mesh_description).at(i).vertices.at(j);

		std::cout << std::endl;
	}

	unsigned int first_cell_shared_line_index = -1;
	unsigned int second_cell_shared_line_index = -1;

	Triangulation<2>::active_cell_iterator cell = triangulation.begin_active();

	for (unsigned int f = 0; f < 4; f++)
		if ( cell->neighbor_index(f) != -1 )
		{
			first_cell_shared_line_index = f;
			break;
		}

	std::cout << "Cell id = 0:" << std::endl;

	std::cout << "The local index of the shared line: " << 
		first_cell_shared_line_index << std::endl;

	std::cout << "Shared line orientation: " << 
		cell->line_orientation(first_cell_shared_line_index) << std::endl;

	cell++;

	std::cout << "Cell id = 1:" << std::endl;

	for (unsigned int f = 0; f < 4; f++)
		if ( cell->neighbor_index(f) != -1 )
		{
			second_cell_shared_line_index = f;
			break;
		}
	std::cout << "The local index of the shared line: " << 
		second_cell_shared_line_index << std::endl;

	std::cout << "Shared line orientation: " <<
		cell->line_orientation(second_cell_shared_line_index) << std::endl;

	std::cout << "------------------------------------- Two-cells mesh\n\n";
}

template<>
void ShapeFunctions<3>::make_mesh()
{
  std::cout << "\n\n";
	std::cout << "         If orientation is 001           \n";
	std::cout << "                                         \n";
	std::cout << "          z                              \n";
 	std::cout << "            ^                            \n";   
 	std::cout << "            |     2                      \n";
 	std::cout << "      6 *---|---*-------* 6 (10)         \n";
 	std::cout << "       /|   |  /|      /|                \n";
 	std::cout << "      / |   | / |  7  / |                \n";
 	std::cout << "   7 /  | 3 |/  |(11)/  |                \n";
	std::cout << "    *   |   *   |   *   |                \n";
 	std::cout << "    |   *---|---*---|---* 5 (9)          \n";
 	std::cout << "    |  / 5  |  / 1  |  /                 \n";
 	std::cout << "    | /     | /     | /                  \n";
 	std::cout << "    |/      |/      |/                   \n";
 	std::cout << "  4 *-------*-------*------>y            \n";
 	std::cout << "           / 0       4 (8)               \n";
 	std::cout << "          /                              \n";
 	std::cout << "         /                               \n";
 	std::cout << "        v x                              \n";
	std::cout << "                                         \n";
	std::cout << "                                         \n";
	std::cout << "                                         \n";
	std::cout << "          z                              \n";
 	std::cout << "            ^                            \n";   
 	std::cout << "            |     2                      \n";
 	std::cout << "      6 *---|---*-------* 6 (10)         \n";
 	std::cout << "       /    |  /       /|                \n";
 	std::cout << "      /     | /       / |                \n";
 	std::cout << "   7 /      |/  7(11)/  |                \n";
 	std::cout << "    *-------*-------*   |                \n";
 	std::cout << "    |     3 |       |   * 5 (9)          \n";
 	std::cout << "    |       |       |  /                 \n";
 	std::cout << "    |       |       | /                  \n";
 	std::cout << "    |       |       |/                   \n";
 	std::cout << "  4 *-------*-------*------>y            \n";
 	std::cout << "           / 0       4 (8)               \n";
 	std::cout << "          /                              \n";
 	std::cout << "         /                               \n";
 	std::cout << "        v x                              \n";
	std::cout << "\n\n";

	std::vector<Point<3>> vertices_left =
	{{ 0.0, 0.0, 0.0 },  // 0
	 {-1.0, 0.0, 0.0 },  // 1
	 {-1.0, 0.0, 1.0 },  // 2
	 { 0.0, 0.0, 1.0 },  // 3
	 { 0.0,-1.0, 0.0 },  // 4
	 {-1.0,-1.0, 0.0 },  // 5
	 {-1.0,-1.0, 1.0 },  // 6
	 { 0.0,-1.0, 1.0 }}; // 7

	std::vector<Point<3>> vertices_right =
	{{ 0.0, 0.0, 0.0 },  // 0
	 {-1.0, 0.0, 0.0 },  // 1
	 {-1.0, 0.0, 1.0 },  // 2
	 { 0.0, 0.0, 1.0 },  // 3
	 { 0.0, 1.0, 0.0 },  // 4
	 {-1.0, 1.0, 0.0 },  // 5
	 {-1.0, 1.0, 1.0 },  // 6
	 { 0.0, 1.0, 1.0 }}; // 7
	
	CellData<3> cell_left, cell_right;

	if (combined_face_orientation % 2 == 0)
	{
		cell_left.vertices = {4, 7, 0, 3, 5, 6, 1, 2};
		
		switch (combined_face_orientation)
		{
			case 0:
				cell_right.vertices = {4, 5, 0, 1, 7, 6, 3, 2}; // 000
				break;
			case 2:
				cell_right.vertices = {7, 4, 3, 0, 6, 5, 2, 1}; // 010 
				break;
			case 4:
				cell_right.vertices = {6, 7, 2, 3, 5, 4, 1, 0}; // 100
				break;
			case 6:
 				cell_right.vertices = {5, 6, 1, 2, 4, 7, 0, 3}; // 110 
				break;
			default:
				break;
		}
	}else
	{
		cell_left.vertices = {4, 0, 5, 1, 7, 3, 6, 2};

		switch (combined_face_orientation)
		{
			case 1:
				cell_right.vertices = {0, 4, 1, 5, 3, 7, 2, 6}; // 001
				break;
			case 3:
				cell_right.vertices = {7, 3, 4, 0, 6, 2, 5, 1}; // 011 
				break;
			case 5:
			 	cell_right.vertices = {6, 2, 7, 3, 5, 1, 4, 0}; // 101
				break;
			case 7:
			 	cell_right.vertices = {5, 1, 6, 2, 4, 0, 7, 3}; // 111
				break;
			default:
				break;
		}
	}

	std::cout << "Two one-cell meshes --------------------------------\n";

	std::cout << "Left-mesh cell  | ";
	for (unsigned j : {0, 1, 2, 3, 4, 5, 6, 7})
		std::cout << std::setw(3) << cell_left.vertices.at(j);
	std::cout << std::endl;

	std::cout << "Right-mesh cell | ";
	for (unsigned j : {0, 1, 2, 3, 4, 5, 6, 7})
		std::cout << std::setw(3) << cell_right.vertices.at(j);
	std::cout << std::endl;
	std::cout << std::endl;

	std::vector<CellData<3>> cell_data_left = {cell_left};
	std::vector<CellData<3>> cell_data_right = {cell_right};

	SubCellData subcell_data;  

	triangulation_left.create_triangulation(vertices_left, cell_data_left, subcell_data);
	triangulation_right.create_triangulation(vertices_right, cell_data_right, subcell_data);
	GridGenerator::merge_triangulations(triangulation_left, triangulation_right, triangulation);

	// Now we have three meshes:
	// triangulation_left - single cell in the third octant; the origin is a vertex. 
	// triangulation_right - single cell in the second octant; the origin is a vertex. 
	// triangulation - exactly the same two cells combined in one mesh. 

	// Save the combined mesh.
	GridOut gridout;
	GridOutFlags::Msh msh_flags(true, true);
	gridout.set_flags(msh_flags);

	std::ofstream ofs("output.msh");
	gridout.write_msh(triangulation, ofs);

  // Let's print diagnostic information on the two-cells mesh, triangulation.
	
	std::cout << "Two-cells mesh -------------------------------------\n";
	
	std::tuple< std::vector<Point<3>>, std::vector<CellData<3>>, SubCellData> mesh_description;
	mesh_description = GridTools::get_coarse_mesh_description(triangulation);

	std::cout << "Vertices:\n";
	for (unsigned int i = 0; i < std::get<0>(mesh_description).size(); i++)
	{
		std::cout << std::setw(4) << i << " | ";

		for (unsigned j : {0, 1, 2})
			std::cout << std::setw(4) << std::right << std::get<0>(mesh_description).at(i)[j];

		std::cout << std::endl;
	}

	std::cout << "Cells:\n";
	for (unsigned int i : {0, 1})
	{
		std::cout << "Cell id = " << i << " | ";
	
		for (unsigned int j = 0; j < 8; j++)
			std::cout << std::setw(3) << std::get<1>(mesh_description).at(i).vertices.at(j);

		std::cout << std::endl;
	}

	unsigned int first_cell_shared_face_index = -1;
	unsigned int second_cell_shared_face_index = -1;

	Triangulation<3>::active_cell_iterator cell = triangulation.begin_active();

	for (unsigned int f = 0; f < 6; f++)
		if ( cell->neighbor_index(f) != -1 )
		{
			first_cell_shared_face_index = f;
			break;
		}

	std::cout << "Cell id = 0:" << std::endl;

	std::cout << "The local index of the shared face: " << 
		first_cell_shared_face_index << std::endl;

	std::cout << "Shared face - flip: " << cell->face_flip(first_cell_shared_face_index) << std::endl;
	std::cout << "Shared face - rotation: " << cell->face_rotation(first_cell_shared_face_index) << std::endl;
	std::cout << "Shared face - orientation: " << cell->face_orientation(first_cell_shared_face_index) << std::endl;
	std::cout << "Shared face - combined orientation: " <<
	static_cast<unsigned int>(internal::combined_face_orientation(
				cell->face_orientation(first_cell_shared_face_index), 
				cell->face_rotation(first_cell_shared_face_index),
				cell->face_flip(first_cell_shared_face_index))) << std::endl;  	

	cell++;

	std::cout << "Cell id = 1:" << std::endl;

	for (unsigned int f = 0; f < 6; f++)
		if ( cell->neighbor_index(f) != -1 )
		{
			second_cell_shared_face_index = f;
			break;
		}
	std::cout << "The local index of the shared face: " << 
		second_cell_shared_face_index << std::endl;

	std::cout << "Shared face - flip: " << cell->face_flip(first_cell_shared_face_index) << std::endl;
	std::cout << "Shared face - rotation: " << cell->face_rotation(first_cell_shared_face_index) << std::endl;
	std::cout << "Shared face - orientation: " << cell->face_orientation(first_cell_shared_face_index) << std::endl;
	std::cout << "Shared face - combined orientation: " <<
	static_cast<unsigned int>(internal::combined_face_orientation(
				cell->face_orientation(first_cell_shared_face_index), 
				cell->face_rotation(first_cell_shared_face_index),
				cell->face_flip(first_cell_shared_face_index))) << std::endl;  	

	std::cout << "------------------------------------- Two-cells mesh\n\n";
}

template<int dim>
void ShapeFunctions<dim>::setup_system_left()
{
	dof_handler_left.reinit(triangulation_left);
	dof_handler_left.distribute_dofs(fe);

	DynamicSparsityPattern dsp(dof_handler_left.n_dofs(), dof_handler_left.n_dofs());
	DoFTools::make_sparsity_pattern(
		dof_handler_left,
		dsp,
		AffineConstraints<double>(),
		false);

	sparsity_pattern_left.copy_from(dsp);
	solution_left.reinit(dof_handler_left.n_dofs());
}

template<int dim>
void ShapeFunctions<dim>::setup_system_right()
{
	dof_handler_right.reinit(triangulation_right);
	dof_handler_right.distribute_dofs(fe);

	DynamicSparsityPattern dsp(dof_handler_right.n_dofs(), dof_handler_right.n_dofs());
	DoFTools::make_sparsity_pattern(
		dof_handler_right,
		dsp,
		AffineConstraints<double>(),
		false);

	sparsity_pattern_right.copy_from(dsp);
	solution_right.reinit(dof_handler_right.n_dofs());
}

template<int dim>
void ShapeFunctions<dim>::setup_system()
{
	dof_handler.reinit(triangulation);
	dof_handler.distribute_dofs(fe);

	DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(
		dof_handler,
		dsp,
		AffineConstraints<double>(),
		false);

	sparsity_pattern.copy_from(dsp);
	solution.reinit(dof_handler.n_dofs());
}

template <int dim>
void ShapeFunctions<dim>::save_left(std::string fname) const
{
	std::vector<std::string> solution_names(dim, "ShapeFunction");
	std::vector<DataComponentInterpretation::DataComponentInterpretation>
	interpretation(dim,
		DataComponentInterpretation::component_is_part_of_vector);

	DataOut<dim> data_out;

	data_out.add_data_vector(dof_handler_left,
		solution_left,
		solution_names,
		interpretation);

	data_out.build_patches(16);

	std::ofstream out(fname);

	data_out.write_vtk(out);
	out.close();
}

template <int dim>
void ShapeFunctions<dim>::save_right(std::string fname) const
{
	std::vector<std::string> solution_names(dim, "ShapeFunction");
	std::vector<DataComponentInterpretation::DataComponentInterpretation>
	interpretation(dim,
		DataComponentInterpretation::component_is_part_of_vector);

	DataOut<dim> data_out;

	data_out.add_data_vector(dof_handler_right,
		solution_right,
		solution_names,
		interpretation);

	data_out.build_patches(16);

	std::ofstream out(fname);

	data_out.write_vtk(out);
	out.close();
}

template <int dim>
void ShapeFunctions<dim>::save(std::string fname) const
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

template<int dim>
void ShapeFunctions<dim>::run()
{
	make_mesh();
	setup_system_left();
	setup_system_right();
	setup_system();

	std::cout << fe.get_name() << std::endl;
	std::cout << "Number of dofs, two-cells mesh: " << dof_handler.n_dofs() << std::endl;
	std::cout << "Number of dofs per cell: " << dof_handler_left.get_fe().n_dofs_per_cell() << std::endl;
	std::cout << "Number of dofs per line: " << dof_handler_left.get_fe().n_dofs_per_line() << std::endl;
	std::cout << "Number of dofs per quad: " << dof_handler_left.get_fe().n_dofs_per_quad() << std::endl;
	std::cout << "Number of dofs per hex: " << dof_handler_left.get_fe().n_dofs_per_hex() << std::endl;

	std::cout << "\nWriting to ./Data/\n...\n";

	unsigned int N = solution_left.size();
	for (unsigned int i = 0; i < N; i++)
	{
		solution_left.reinit(N);
		solution_left(i) = 1.0;
		save_left(fname_vtk +"_left_" + std::to_string(i) + ".vtk");

		solution_right.reinit(N);
		solution_right(i) = 1.0;
		save_right(fname_vtk +"_right_" + std::to_string(i) + ".vtk");
	}

	N = solution.size();
	for (unsigned int i = 0; i < N; i++)
	{
		solution.reinit(N);
		solution(i) = 1.0;
		save(fname_vtk +"_both_" + std::to_string(i) + ".vtk");
	}
}

int main()
{
	ShapeFunctions<DIMENSION__> shape_functions;

	return 0;
}

