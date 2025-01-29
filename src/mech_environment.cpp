#include "mech_environment.h"

using namespace biofvm;
using namespace micromech;

mech_environment::mech_environment(microenvironment& m, real_t timestep, real_t agent_types_count)
	: m(m), timestep(timestep), agent_types_count(agent_types_count), agent_data(*this)
{}
