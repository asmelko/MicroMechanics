#include "base_wall_membrane_model.h"

#include <base_membrane_data.h>

#include "mech_environment.h"

using namespace micromech;
using namespace biofvm;

constexpr void update_membrane_velocity(real_t position, real_t bounding_box, real_t sign, real_t radius,
										real_t repulsion_strength, real_t& velocity)
{
	biofvm::real_t distance = std::abs(bounding_box - position);

	distance = std::max<biofvm::real_t>(distance, 0.00001);

	biofvm::real_t repulsion = 1 - distance / radius;
	repulsion = repulsion < 0 ? 0 : repulsion;

	repulsion *= repulsion * repulsion_strength * sign;

	velocity += repulsion * distance;
}

template <index_t dims>
constexpr void update_membrane_velocities(biofvm::real_t* __restrict__ velocity,
										  const biofvm::real_t* __restrict__ position,
										  const biofvm::cartesian_mesh& mesh, const biofvm::real_t radius,
										  const biofvm::real_t repulsion_strength)
{
	if constexpr (dims == 1)
	{
		update_membrane_velocity(position[0], mesh.bounding_box_mins[0], 1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[0], mesh.bounding_box_maxs[1], -1, radius, repulsion_strength, velocity[0]);
	}
	if constexpr (dims == 2)
	{
		update_membrane_velocity(position[0], mesh.bounding_box_mins[0], 1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[0], mesh.bounding_box_maxs[0], -1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[1], mesh.bounding_box_mins[1], 1, radius, repulsion_strength, velocity[1]);
		update_membrane_velocity(position[1], mesh.bounding_box_maxs[1], -1, radius, repulsion_strength, velocity[1]);
	}
	if constexpr (dims <= 3)
	{
		update_membrane_velocity(position[0], mesh.bounding_box_mins[0], 1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[0], mesh.bounding_box_maxs[0], -1, radius, repulsion_strength, velocity[0]);
		update_membrane_velocity(position[1], mesh.bounding_box_mins[1], 1, radius, repulsion_strength, velocity[1]);
		update_membrane_velocity(position[1], mesh.bounding_box_maxs[1], -1, radius, repulsion_strength, velocity[1]);
		update_membrane_velocity(position[2], mesh.bounding_box_mins[2], 1, radius, repulsion_strength, velocity[2]);
		update_membrane_velocity(position[2], mesh.bounding_box_maxs[2], -1, radius, repulsion_strength, velocity[2]);
	}
}

template <index_t dims>
void update_basement_membrane_interactions_internal(index_t agents_count, real_t* __restrict__ velocity,
													const real_t* __restrict__ position,
													const real_t* __restrict__ radius,
													const real_t* __restrict__ cell_BM_repulsion_strength,
													const std::uint8_t* __restrict__ is_movable,
													const cartesian_mesh& mesh)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		if (is_movable[i] == 0)
			continue;

		update_membrane_velocities<dims>(velocity + i * dims, position + i * dims, mesh, radius[i],
										 cell_BM_repulsion_strength[i]);
	}
}

void base_wall_membrane_model::compute_basement_membrane_interactions(mech_environment& me)
{
	auto& data = me.agent_data;
	auto& membrane_data = static_cast<base_membrane_data&>(*data.membrane_data.get());

	if (me.m.mesh.dims == 1)
		update_basement_membrane_interactions_internal<1>(
			data.bio_agent_data.agents_count, data.velocity.data(), data.bio_agent_data.positions.data(),
			data.radius.data(), membrane_data.cell_BM_repulsion_strength.data(), data.is_movable.data(), me.m.mesh);
	else if (me.m.mesh.dims == 2)
		update_basement_membrane_interactions_internal<2>(
			data.bio_agent_data.agents_count, data.velocity.data(), data.bio_agent_data.positions.data(),
			data.radius.data(), membrane_data.cell_BM_repulsion_strength.data(), data.is_movable.data(), me.m.mesh);
	else if (me.m.mesh.dims == 3)
		update_basement_membrane_interactions_internal<3>(
			data.bio_agent_data.agents_count, data.velocity.data(), data.bio_agent_data.positions.data(),
			data.radius.data(), membrane_data.cell_BM_repulsion_strength.data(), data.is_movable.data(), me.m.mesh);
}

base_wall_membrane_model::base_wall_membrane_model(mech_environment& me)
{
	if (dynamic_cast<base_membrane_data*>(me.agent_data.membrane_data.get()) == nullptr)
	{
		me.agent_data.membrane_data = std::make_unique<base_membrane_data>(me);
	}
}
