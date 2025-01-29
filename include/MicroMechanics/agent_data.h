#pragma once

#include <BioFVM/types.h>

namespace micromech {

struct agent_data
{
	virtual void add(biofvm::index_t size, biofvm::index_t cell_definitions_count);
	virtual void remove(biofvm::index_t index, biofvm::index_t size, biofvm::index_t cell_definitions_count);
};

} // namespace micromech
