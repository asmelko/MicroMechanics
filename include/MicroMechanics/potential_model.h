#pragma once

#include "mech_agent_data.h"

namespace micromech {

class potential_model
{
public:
	virtual void update_velocities(mech_agent_data& cells) = 0;

	virtual void update_neighbors(mech_agent_data& cells) = 0;
};

} // namespace micromech
