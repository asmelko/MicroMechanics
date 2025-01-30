#pragma once

#include <vector>

#include "agent_data.h"

namespace micromech {

struct base_membrane_data : public agent_data
{
	// std::vector<biofvm::real_t> cell_BM_adhesion_strength;
	std::vector<biofvm::real_t> cell_BM_repulsion_strength;

	base_membrane_data(mech_environment& me);

	virtual void add() override;
	virtual void remove(biofvm::index_t index) override;
};

} // namespace micromech
