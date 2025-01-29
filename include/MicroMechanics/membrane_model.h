#pragma once

#include "mech_environment.h"

namespace micromech {

class membrane_model
{
public:
	virtual void compute_basement_membrane_interactions(mech_environment& me) = 0;
};

} // namespace micromech
