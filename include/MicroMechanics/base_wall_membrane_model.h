#pragma once

#include "mech_environment.h"
#include "membrane_model.h"

namespace micromech {

class base_wall_membrane_model : public membrane_model
{
public:
	base_wall_membrane_model(mech_environment& me);
	
	virtual void compute_basement_membrane_interactions(mech_environment& me) = 0;
};

} // namespace micromech
