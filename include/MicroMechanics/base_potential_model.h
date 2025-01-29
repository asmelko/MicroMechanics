#pragma once

#include <BioFVM/mesh.h>
#include <BioFVM/types.h>

#include "grid_space_partitioner.h"
#include "mech_environment.h"
#include "potential_model.h"

namespace micromech {

class base_potential_model : public potential_model
{
	void compute_agents_potentials(mech_environment& me);
	void attach_detach_springs(mech_environment& me);
	void compute_springs_potentials(mech_environment& me);

	grid_space_partitioner& partitioner_;

public:
	base_potential_model(grid_space_partitioner& partitioner, mech_environment& me);

	virtual void update_velocities(mech_environment& me) override;

	virtual void update_neighbors(mech_environment& me) override;

	virtual void update_positions(mech_environment& me) override;
};

} // namespace micromech
