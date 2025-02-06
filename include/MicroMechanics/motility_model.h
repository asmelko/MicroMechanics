#pragma once

namespace micromech {

struct mech_environment;

class motility_model
{
public:
	virtual void update_motility_velocities(mech_environment& me) = 0;
};

} // namespace micromech
