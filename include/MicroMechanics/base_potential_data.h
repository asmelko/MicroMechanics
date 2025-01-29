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

	virtual void add(biofvm::index_t size, biofvm::index_t cell_definitions_count) override;
	virtual void remove(biofvm::index_t index, biofvm::index_t size, biofvm::index_t cell_definitions_count) override;
};

} // namespace micromech
