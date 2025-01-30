#pragma once

#include <cstdint>
#include <memory>

#include <BioFVM/agent_data.h>

#include "agent_data.h"

namespace micromech {

struct mech_environment;

struct mech_agent_data
{
	biofvm::agent_data bio_agent_data;

	mech_environment& me;

	std::vector<biofvm::real_t> velocity;
	std::vector<biofvm::real_t> radius;
	std::vector<std::uint8_t> is_movable;

	std::vector<biofvm::index_t> agent_type_indices;

	std::vector<std::vector<biofvm::index_t>> neighbors;

	std::unique_ptr<agent_data> potential_data, membrane_data, motility_data;

	mech_agent_data(mech_environment& me);

	void add();
	void remove(biofvm::index_t index);

	biofvm::index_t agents_count() const;
};

} // namespace micromech
