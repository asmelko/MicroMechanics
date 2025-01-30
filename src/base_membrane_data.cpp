#include "base_membrane_data.h"

using namespace biofvm;
using namespace micromech;

base_membrane_data::base_membrane_data(mech_environment& me) : agent_data(me) {}

void base_membrane_data::add() { cell_BM_repulsion_strength.resize(agents_count(), 0); }

void base_membrane_data::remove(index_t index)
{
	if (index == agents_count())
		return;

	cell_BM_repulsion_strength[index] = cell_BM_repulsion_strength[agents_count()];
}
