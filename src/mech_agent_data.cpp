#include "mech_agent_data.h"

#include <memory>

#include <BioFVM/data_utils.h>

#include "empty_data.h"
#include "mech_environment.h"

using namespace biofvm;
using namespace micromech;

mech_agent_data::mech_agent_data(mech_environment& me)
	: bio_agent_data(me.m),
	  me(me),
	  potential_data(std::make_unique<empty_data>(me)),
	  membrane_data(std::make_unique<empty_data>(me)),
	  motility_data(std::make_unique<empty_data>(me))
{}

void mech_agent_data::add()
{
	bio_agent_data.add();
	potential_data->add();
	membrane_data->add();
	motility_data->add();

	velocity.resize(agents_count() * me.m.mesh.dims);
	radius.resize(agents_count());
	is_movable.resize(agents_count());
	agent_type_indices.resize(agents_count());
	neighbors.resize(agents_count());
}

void mech_agent_data::remove(index_t index)
{
	bio_agent_data.remove(index);
	potential_data->remove(index);
	membrane_data->remove(index);
	motility_data->remove(index);

	if (index == agents_count())
		return;

	move_vector(velocity.data() + index * me.m.mesh.dims, velocity.data() + agents_count() * me.m.mesh.dims,
				me.m.mesh.dims);

	radius[index] = radius[agents_count()];
	is_movable[index] = is_movable[agents_count()];
	agent_type_indices[index] = agent_type_indices[agents_count()];
	neighbors[index] = neighbors[agents_count()];
}
