#include "agent_data.h"

#include "mech_environment.h"

micromech::agent_data::agent_data(micromech::mech_environment& me) : me(me) {}

biofvm::index_t micromech::agent_data::agents_count() const { return me.agent_data.agents_count(); }
