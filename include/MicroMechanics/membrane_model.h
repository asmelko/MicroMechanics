#pragma once

namespace micromech {

struct mech_environment;

class membrane_model
{
public:
	virtual void compute_basement_membrane_interactions(mech_environment& me) = 0;
};

} // namespace micromech
