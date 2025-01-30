#include "base_motility_data.h"

#include <BioFVM/data_utils.h>

#include "mech_environment.h"

using namespace biofvm;
using namespace micromech;

base_motility_data::base_motility_data(mech_environment& me) : agent_data(me) {}

void base_motility_data::add()
{
	is_motile.resize(agents_count());
	persistence_time.resize(agents_count());
	migration_speed.resize(agents_count());

	migration_bias_direction.resize(agents_count() * me.m.mesh.dims);
	migration_bias.resize(agents_count());

	motility_vector.resize(agents_count() * me.m.mesh.dims);

	restrict_to_2d.resize(agents_count());

	chemotaxis_index.resize(agents_count());
	chemotaxis_direction.resize(agents_count());
	chemotactic_sensitivities.resize(agents_count() * me.m.substrates_count);

	update_migration_bias_direction.resize(agents_count());
}

void base_motility_data::remove(index_t index)
{
	if (index == agents_count())
		return;
    
	is_motile[index] = is_motile[agents_count()];
	persistence_time[index] = persistence_time[agents_count()];
	migration_speed[index] = migration_speed[agents_count()];

	move_vector(migration_bias_direction.data() + index * me.m.mesh.dims,
				migration_bias_direction.data() + agents_count() * me.m.mesh.dims, me.m.mesh.dims);
	migration_bias[index] = migration_bias[agents_count()];

	move_vector(motility_vector.data() + index * me.m.mesh.dims,
				motility_vector.data() + agents_count() * me.m.mesh.dims, me.m.mesh.dims);

	restrict_to_2d[index] = restrict_to_2d[agents_count()];

	chemotaxis_index[index] = chemotaxis_index[agents_count()];
	chemotaxis_direction[index] = chemotaxis_direction[agents_count()];
	move_vector(chemotactic_sensitivities.data() + index * me.m.substrates_count,
				chemotactic_sensitivities.data() + agents_count() * me.m.substrates_count, me.m.substrates_count);

	update_migration_bias_direction[index] = update_migration_bias_direction[agents_count()];
}
