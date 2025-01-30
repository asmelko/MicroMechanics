#pragma once

#include <BioFVM/types.h>

namespace micromech {

struct mech_environment;

struct agent_data
{
	mech_environment& me;

	agent_data(mech_environment& me);

	biofvm::index_t agents_count() const;

	virtual void add() = 0;
	virtual void remove(biofvm::index_t index) = 0;
};

} // namespace micromech
