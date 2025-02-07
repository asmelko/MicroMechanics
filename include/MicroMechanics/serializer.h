#pragma once

#include "mech_environment.h"

namespace micromech {

class serializer
{
public:
	virtual void serialize_one_timestep(const mech_environment& me) = 0;
};

} // namespace micromech