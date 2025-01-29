#pragma once

#include "mech_agent_data.h"

namespace micromech {

class membrane_model
{
public:
	virtual void compute_basement_membrane_interactions(mech_agent_data& cells) = 0;
};

} // namespace micromech
