#pragma once

#include <cstdint>
#include <functional>

#include <BioFVM/agent_data.h>

#include "agent_data.h"

namespace micromech {

struct base_motility_data : public agent_data
{
	using direction_update_func = std::function<void(biofvm::real_t*)>;

	std::vector<std::uint8_t> is_motile;
	std::vector<biofvm::real_t> persistence_time;
	std::vector<biofvm::real_t> migration_speed;

	std::vector<biofvm::real_t> migration_bias_direction;
	std::vector<biofvm::real_t> migration_bias;

	std::vector<biofvm::real_t> motility_vector;

	std::vector<std::uint8_t> restrict_to_2d;

	std::vector<biofvm::index_t> chemotaxis_index;
	std::vector<biofvm::index_t> chemotaxis_direction;
	std::vector<biofvm::real_t> chemotactic_sensitivities;

	std::vector<direction_update_func> update_migration_bias_direction;

	base_motility_data(mech_environment& me);

	virtual void add() override;
	virtual void remove(biofvm::index_t index) override;
};

} // namespace micromech
