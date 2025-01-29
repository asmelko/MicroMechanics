#pragma once

#include "membrane_model.h"

namespace micromech {

class base_wall_membrane_model : public membrane_model
{
public:
	virtual void compute_basement_membrane_interactions(mech_agent_data& cells) = 0;
};

} // namespace micromech
