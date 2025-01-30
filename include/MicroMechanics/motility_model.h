#pragma once

#include "mech_environment.h"

namespace micromech {

class motility_model
{
public:
	virtual void update_motility_velocities(mech_environment& me) = 0;
};

} // namespace micromech
