#pragma once

namespace micromech {

struct mech_environment;

class potential_model
{
public:
	virtual void update_velocities(mech_environment& me) = 0;

	virtual void update_neighbors(mech_environment& me) = 0;

	virtual void update_positions(mech_environment& me) = 0;
};

} // namespace micromech
