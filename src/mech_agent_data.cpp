#include "mech_agent_data.h"

#include <memory>

#include "empty_data.h"
#include "mech_environment.h"

using namespace biofvm;
using namespace micromech;

mech_agent_data::mech_agent_data(mech_environment& me)
	: bio_agent_data(me.m),
	  me(me),
	  potential_data(std::make_unique<empty_data>(me)),
	  membrane_data(std::make_unique<empty_data>(me)),
	  motility_data(std::make_unique<empty_data>(me))
{}
