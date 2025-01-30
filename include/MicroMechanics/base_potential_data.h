#pragma once

#include <BioFVM/agent_data.h>

#include "agent_data.h"

namespace micromech {

struct base_potential_data : public agent_data
{
	std::vector<biofvm::real_t> cell_cell_adhesion_strength;
	std::vector<biofvm::real_t> cell_cell_repulsion_strength;

	std::vector<biofvm::real_t> cell_adhesion_affinities;

	std::vector<biofvm::real_t> relative_maximum_adhesion_distance;

	std::vector<biofvm::index_t> maximum_number_of_attachments;
	std::vector<biofvm::real_t> attachment_elastic_constant;

	std::vector<biofvm::real_t> attachment_rate;
	std::vector<biofvm::real_t> detachment_rate;

	std::vector<biofvm::real_t> simple_pressure;

	std::vector<biofvm::real_t> previous_velocity;

	std::vector<std::vector<biofvm::index_t>> springs;

	base_potential_data(mech_environment& me);

	virtual void add() override;
	virtual void remove(biofvm::index_t index) override;
};

} // namespace micromech
