#pragma once

#include <cstdint>
#include <memory>

#include <BioFVM/agent_data.h>

#include "agent_data.h"

namespace micromech {

struct mech_agent_data
{
	biofvm::agent_data bio_agent_data;

	std::vector<biofvm::real_t> velocities;
	std::vector<biofvm::real_t> radii;
	std::vector<std::uint8_t> movable;

	std::unique_ptr<agent_data> potential_data, membrane_data;
};

} // namespace micromech
