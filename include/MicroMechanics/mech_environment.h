#pragma once

#include <memory>

#include <BioFVM/microenvironment.h>

#include "mech_agent_data.h"
#include "membrane_model.h"
#include "motility_model.h"
#include "potential_model.h"

namespace micromech {

struct mech_environment
{
	biofvm::microenvironment& m;

	biofvm::real_t timestep;

	biofvm::index_t agent_types_count;

	mech_agent_data agent_data;

	std::unique_ptr<potential_model> potential_m;
	std::unique_ptr<membrane_model> membrane_m;
	std::unique_ptr<motility_model> motility_m;

	mech_environment(biofvm::microenvironment& m, biofvm::real_t timestep, biofvm::real_t agent_types_count);
};

} // namespace micromech
