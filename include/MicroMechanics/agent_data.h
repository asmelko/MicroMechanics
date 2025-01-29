#pragma once

#include <BioFVM/types.h>

namespace micromech {

struct mech_environment;

struct agent_data
{
	mech_environment& me;

	agent_data(mech_environment& me);

	virtual void add();
	virtual void remove(biofvm::index_t index);
};

} // namespace micromech
