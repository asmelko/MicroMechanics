#pragma once

#include "agent_data.h"

namespace micromech {

struct empty_data : public agent_data
{
	empty_data(mech_environment& me) : agent_data(me) {}

	virtual void add() override {}
	virtual void remove(biofvm::index_t) override {}
};

} // namespace micromech
