#pragma once

#include <memory>
#include <vector>

#include <BioFVM/mesh.h>
#include <BioFVM/types.h>

#include "potential_model.h"

namespace micromech {

class base_potential_model : public potential_model
{
	void compute_agents_potentials();
	void attach_detach_springs();
	void compute_springs_potentials();

protected:
	biofvm::real_t timestep;
	biofvm::cartesian_mesh mechanics_mesh;
	std::unique_ptr<std::vector<biofvm::index_t>> agents_in_voxels;

public:
	base_potential_model(const biofvm::real_t timestep);

	virtual void update_velocities(mech_agent_data& cells) override;

	virtual void update_neighbors(mech_agent_data& cells) override;
};

} // namespace micromech
