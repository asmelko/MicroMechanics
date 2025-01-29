#pragma once

#include "mech_environment.h"
#include "membrane_model.h"

namespace micromech {

class no_wall_membrane_model : public membrane_model
{
public:
	virtual void compute_basement_membrane_interactions(mech_environment&) {}
};

} // namespace micromech
