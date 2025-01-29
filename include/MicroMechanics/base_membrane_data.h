#pragma once
#include <vector>

#include "agent_data.h"

namespace micromech {

struct base_membrane_data : public agent_data
{
	std::vector<biofvm::real_t> cell_BM_adhesion_strength;
	std::vector<biofvm::real_t> cell_BM_repulsion_strength;

	virtual void add(biofvm::index_t size, biofvm::index_t cell_definitions_count) override;
	virtual void remove(biofvm::index_t index, biofvm::index_t size, biofvm::index_t cell_definitions_count) override;
};

} // namespace micromech
