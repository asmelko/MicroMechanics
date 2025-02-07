#pragma once

#include <string_view>

#include <BioFVM/vtk_serializer_base.h>


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#pragma GCC diagnostic pop

#include "base_potential_data.h"
#include "serializer.h"


namespace micromech {

class base_vtk_serializer : public serializer, public biofvm::vtk_serializer_base
{
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer_;

	vtkSmartPointer<vtkPoints> points_;
	vtkSmartPointer<vtkFloatArray> point_data_;
	std::vector<vtkSmartPointer<vtkLine>> lines_;


	void serialize_potential_data(const base_potential_data& data, vtkSmartPointer<vtkUnstructuredGrid>& grid);

protected:
	void serialize_cell_data(const mech_environment& me, vtkSmartPointer<vtkUnstructuredGrid>& grid);

public:
	base_vtk_serializer(std::string_view output_dir);

	virtual void serialize_one_timestep(const mech_environment& m) override;
};

} // namespace micromech
