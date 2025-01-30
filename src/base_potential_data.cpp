#include "base_potential_data.h"

#include <BioFVM/data_utils.h>

#include "mech_environment.h"

using namespace biofvm;
using namespace micromech;

base_potential_data::base_potential_data(mech_environment& me) : agent_data(me) {}

void base_potential_data::add()
{
	cell_cell_adhesion_strength.resize(agents_count());
	cell_cell_repulsion_strength.resize(agents_count());

	cell_adhesion_affinities.resize(agents_count() * me.agent_types_count);

	relative_maximum_adhesion_distance.resize(agents_count());

	maximum_number_of_attachments.resize(agents_count());
	attachment_elastic_constant.resize(agents_count());

	attachment_rate.resize(agents_count());
	detachment_rate.resize(agents_count());

	simple_pressure.resize(agents_count());

	previous_velocity.resize(agents_count() * me.m.mesh.dims);

	springs.resize(agents_count());
}

void base_potential_data::remove(index_t index)
{
	if (index == agents_count())
		return;

	cell_cell_adhesion_strength[index] = cell_cell_adhesion_strength[agents_count()];
	cell_cell_repulsion_strength[index] = cell_cell_repulsion_strength[agents_count()];

	move_vector(cell_adhesion_affinities.data() + index * me.agent_types_count,
				cell_adhesion_affinities.data() + agents_count() * me.agent_types_count, me.agent_types_count);

	relative_maximum_adhesion_distance[index] = relative_maximum_adhesion_distance[agents_count()];

	maximum_number_of_attachments[index] = maximum_number_of_attachments[agents_count()];
	attachment_elastic_constant[index] = attachment_elastic_constant[agents_count()];

	attachment_rate[index] = attachment_rate[agents_count()];
	detachment_rate[index] = detachment_rate[agents_count()];

	simple_pressure[index] = simple_pressure[agents_count()];

	move_vector(previous_velocity.data() + index * me.m.mesh.dims,
				previous_velocity.data() + agents_count() * me.m.mesh.dims, me.m.mesh.dims);

	springs[index] = springs[agents_count()];
}
