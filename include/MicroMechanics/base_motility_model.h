#pragma once

#include <BioFVM/mesh.h>
#include <BioFVM/types.h>

#include "mech_environment.h"
#include "motility_model.h"

namespace micromech {

class base_motility_model : public motility_model
{
public:
	base_motility_model(mech_environment& me);

	virtual void update_motility_velocities(mech_environment& me) override;
};

} // namespace micromech
