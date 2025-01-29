#pragma once

#include <BioFVM/microenvironment.h>

#include "mech_agent_data.h"

namespace micromech {

struct mech_environment
{
	biofvm::microenvironment& m;

	biofvm::real_t timestep;

	biofvm::real_t agent_types_count;

    mech_agent_data agent_data;

	mech_environment(biofvm::microenvironment& m, biofvm::real_t timestep, biofvm::real_t agent_types_count);
};

} // namespace micromech
