#include "base_vtk_serializer.h"

#include <filesystem>

#include "base_potential_data.h"
#include "mech_agent_data.h"
#include "mech_environment.h"

using namespace biofvm;
using namespace micromech;

void base_vtk_serializer::serialize_potential_data(const base_potential_data&, vtkSmartPointer<vtkUnstructuredGrid>&) {}

point_t<real_t, 3> get_tuple(index_t i, index_t dims, const real_t* data)
{
	if (dims == 1)
	{
		return point_t<real_t, 3> { data[i], 0, 0 };
	}
	if (dims == 2)
	{
		return point_t<real_t, 3> { data[i * 2], data[i * 2 + 1], 0 };
	}
	else
	{
		return point_t<real_t, 3> { data[i * 3], data[i * 3 + 1], data[i * 3 + 2] };
	}
}

void base_vtk_serializer::serialize_cell_data(const mech_environment& me, vtkSmartPointer<vtkUnstructuredGrid>& grid)
{
	points_->SetNumberOfPoints(me.agent_data.agents_count());
#ifdef USE_DOUBLES
	points->SetDataTypeToDouble();
#else
	points_->SetDataTypeToFloat();
#endif

	for (index_t cell_idx = 0; cell_idx < me.agent_data.agents_count(); ++cell_idx)
	{
		point_t<real_t, 3> pos = get_tuple(cell_idx, me.m.mesh.dims, me.agent_data.bio_agent_data.positions.data());
		points_->InsertPoint(cell_idx, pos.data());
	}

	grid->SetPoints(points_);

	point_data_->SetNumberOfComponents(3);
	point_data_->SetNumberOfTuples(me.agent_data.agents_count());
	point_data_->SetName("velocity");

	auto& velocity = dynamic_cast<base_potential_data*>(me.agent_data.potential_data.get())->previous_velocity;
	for (index_t cell_idx = 0; cell_idx < me.agent_data.agents_count(); ++cell_idx)
	{
		point_t<real_t, 3> vel = get_tuple(cell_idx, me.m.mesh.dims, velocity.data());
		point_data_->SetTuple(cell_idx, vel.data());
	}

	grid->GetPointData()->AddArray(point_data_);

	lines_.resize(0);

	for (index_t cell_idx = 0; cell_idx < me.agent_data.agents_count(); ++cell_idx)
	{
		for (index_t nei_idx : me.agent_data.neighbors[cell_idx])
		{
			if (nei_idx > cell_idx)
				continue;

			lines_.push_back(vtkSmartPointer<vtkLine>::New());
			lines_.back()->GetPointIds()->SetId(0, cell_idx);
			lines_.back()->GetPointIds()->SetId(1, nei_idx);
			grid->InsertNextCell(lines_.back()->GetCellType(), lines_.back()->GetPointIds());
		}
	}
}

base_vtk_serializer::base_vtk_serializer(std::string_view output_dir)
	: vtk_serializer_base(output_dir, "vtk_micromechanics"),
	  writer_(vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New()),
	  points_(vtkSmartPointer<vtkPoints>::New()),
	  point_data_(vtkSmartPointer<vtkFloatArray>::New())
{
	// writer_->SetDataModeToBinary();
	// writer_->SetCompressorTypeToLZ4();

	writer_->SetDataModeToAscii();
	writer_->SetCompressorTypeToNone();
}

void base_vtk_serializer::serialize_one_timestep(const mech_environment& me)
{
	vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	serialize_cell_data(me, grid);

	std::ostringstream ss;

	ss << "micromechanics_" << std::setw(6) << std::setfill('0') << iteration_ << ".vtu";

	auto file_name = ss.str();
	auto file_path = std::filesystem::path(vtks_dir_) / file_name;

	writer_->SetInputData(grid);
	writer_->SetFileName(file_path.string().c_str());
	writer_->Write();

	append_to_pvd(file_name);
	++iteration_;
}
