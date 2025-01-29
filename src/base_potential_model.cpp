#include "base_potential_model.h"

#include <BioFVM/microenvironment.h>

#include "base_potential_data.h"
#include "grid_space_partitioner.h"
#include "mech_environment.h"
#include "potentials_helper.h"
#include "random.h"

using namespace micromech;
using namespace biofvm;

base_potential_model::base_potential_model(grid_space_partitioner& partitioner, mech_environment& me)
	: partitioner_(partitioner)
{
	if (dynamic_cast<base_potential_data*>(me.agent_data.potential_data.get()) == nullptr)
	{
		me.agent_data.potential_data = std::make_unique<base_potential_data>(me);
	}
}

template <index_t dims>
void update_cell_neighbors_internal(index_t agents_count, const real_t* __restrict__ position,
									const real_t* __restrict__ radius,
									const real_t* __restrict__ relative_maximum_adhesion_distance,
									const std::uint8_t* __restrict__ is_movable,
									std::vector<index_t>* __restrict__ neighbors, grid_space_partitioner& partitioner)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		if (is_movable[i] == 0)
			continue;

		partitioner.for_each_in_neighborhood<dims>(position + dims * i, i, [=](index_t j) {
			const real_t adhesion_distance =
				relative_maximum_adhesion_distance[i] * radius[i] + relative_maximum_adhesion_distance[j] * radius[j];

			const real_t distance = potentials_helper<dims>::distance(position + i * dims, position + j * dims);

			if (distance <= adhesion_distance)
			{
				neighbors[i].push_back(j);
			}
		});
	}
}

void base_potential_model::update_neighbors(mech_environment& me)
{
	auto& data = me.agent_data;
	auto& potential_data = static_cast<base_potential_data&>(*data.potential_data.get());

	// clear neighbors
#pragma omp for
	for (index_t i = 0; i < data.agents_count(); i++)
		data.neighbors[i].clear();

	if (me.m.mesh.dims == 1)
		update_cell_neighbors_internal<1>(data.agents_count(), data.bio_agent_data.positions.data(), data.radius.data(),
										  potential_data.relative_maximum_adhesion_distance.data(),
										  data.is_movable.data(), data.neighbors.data(), partitioner_);
	else if (me.m.mesh.dims == 2)
		update_cell_neighbors_internal<2>(data.agents_count(), data.bio_agent_data.positions.data(), data.radius.data(),
										  potential_data.relative_maximum_adhesion_distance.data(),
										  data.is_movable.data(), data.neighbors.data(), partitioner_);
	else if (me.m.mesh.dims == 3)
		update_cell_neighbors_internal<3>(data.agents_count(), data.bio_agent_data.positions.data(), data.radius.data(),
										  potential_data.relative_maximum_adhesion_distance.data(),
										  data.is_movable.data(), data.neighbors.data(), partitioner_);
}

void clear_simple_pressure(real_t* __restrict__ simple_pressure, index_t count)
{
#pragma omp for
	for (index_t i = 0; i < count; i++)
	{
		simple_pressure[i] = 0;
	}
}

template <index_t dims>
void solve_pair(index_t lhs, index_t rhs, index_t cell_defs_count, real_t* __restrict__ velocity,
				real_t* __restrict__ simple_pressure, const real_t* __restrict__ position,
				const real_t* __restrict__ radius, const real_t* __restrict__ cell_cell_repulsion_strength,
				const real_t* __restrict__ cell_cell_adhesion_strength,
				const real_t* __restrict__ relative_maximum_adhesion_distance,
				const real_t* __restrict__ cell_adhesion_affinity, const index_t* __restrict__ cell_definition_index)
{
	constexpr real_t simple_pressure_coefficient = 36.64504274775163; // 1 / (12 * (1 - sqrt(pi/(2*sqrt(3))))^2)

	real_t position_difference[dims];

	const real_t distance = std::max<real_t>(potentials_helper<dims>::difference_and_distance(
												 position + lhs * dims, position + rhs * dims, position_difference),
											 0.00001);

	// compute repulsion
	real_t repulsion;
	{
		const real_t repulsive_distance = radius[lhs] + radius[rhs];

		repulsion = 1 - distance / repulsive_distance;

		repulsion = repulsion < 0 ? 0 : repulsion;

		repulsion *= repulsion;

		// update simple pressure
		simple_pressure[lhs] += repulsion * simple_pressure_coefficient;
		simple_pressure[rhs] += repulsion * simple_pressure_coefficient;

		repulsion *= std::sqrt(cell_cell_repulsion_strength[lhs] * cell_cell_repulsion_strength[rhs]);
	}

	// compute adhesion
	real_t adhesion;
	{
		const real_t adhesion_distance = relative_maximum_adhesion_distance[lhs] * radius[lhs]
										 + relative_maximum_adhesion_distance[rhs] * radius[rhs];

		adhesion = 1 - distance / adhesion_distance;

		adhesion *= adhesion;

		const index_t lhs_cell_def_index = cell_definition_index[lhs];
		const index_t rhs_cell_def_index = cell_definition_index[rhs];

		adhesion *= std::sqrt(cell_cell_adhesion_strength[lhs] * cell_cell_adhesion_strength[rhs]
							  * cell_adhesion_affinity[lhs * cell_defs_count + rhs_cell_def_index]
							  * cell_adhesion_affinity[rhs * cell_defs_count + lhs_cell_def_index]);
	}

	real_t force = (repulsion - adhesion) / distance;

	potentials_helper<dims>::update_velocity(velocity + lhs * dims, position_difference, force);
}

template <index_t dims>
void update_cell_forces_internal(
	index_t agents_count, index_t cell_def_count, real_t* __restrict__ velocity, real_t* __restrict__ simple_pressure,
	const real_t* __restrict__ position, const real_t* __restrict__ radius,
	const real_t* __restrict__ cell_cell_repulsion_strength, const real_t* __restrict__ cell_cell_adhesion_strength,
	const real_t* __restrict__ relative_maximum_adhesion_distance, const index_t* __restrict__ cell_definition_index,
	const real_t* __restrict__ cell_adhesion_affinities, const std::uint8_t* __restrict__ is_movable,
	std::vector<index_t>* __restrict__ neighbors)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		if (is_movable[i] == 0)
			continue;

		for (const index_t j : neighbors[i])
		{
			solve_pair<dims>(i, j, cell_def_count, velocity, simple_pressure, position, radius,
							 cell_cell_repulsion_strength, cell_cell_adhesion_strength,
							 relative_maximum_adhesion_distance, cell_adhesion_affinities, cell_definition_index);
		}
	}
}

void base_potential_model::compute_agents_potentials(mech_environment& me)
{
	auto& data = me.agent_data;
	auto& potential_data = static_cast<base_potential_data&>(*data.potential_data.get());

	clear_simple_pressure(potential_data.simple_pressure.data(), data.agents_count());

	if (me.m.mesh.dims == 1)
		update_cell_forces_internal<1>(
			data.agents_count(), me.agent_types_count, data.velocity.data(), potential_data.simple_pressure.data(),
			data.bio_agent_data.positions.data(), data.radius.data(),
			potential_data.cell_cell_repulsion_strength.data(), potential_data.cell_cell_adhesion_strength.data(),
			potential_data.relative_maximum_adhesion_distance.data(), data.agent_type_indices.data(),
			potential_data.cell_adhesion_affinities.data(), data.is_movable.data(), data.neighbors.data());
	else if (me.m.mesh.dims == 2)
		update_cell_forces_internal<2>(
			data.agents_count(), me.agent_types_count, data.velocity.data(), potential_data.simple_pressure.data(),
			data.bio_agent_data.positions.data(), data.radius.data(),
			potential_data.cell_cell_repulsion_strength.data(), potential_data.cell_cell_adhesion_strength.data(),
			potential_data.relative_maximum_adhesion_distance.data(), data.agent_type_indices.data(),
			potential_data.cell_adhesion_affinities.data(), data.is_movable.data(), data.neighbors.data());
	else if (me.m.mesh.dims == 3)
		update_cell_forces_internal<3>(
			data.agents_count(), me.agent_types_count, data.velocity.data(), potential_data.simple_pressure.data(),
			data.bio_agent_data.positions.data(), data.radius.data(),
			potential_data.cell_cell_repulsion_strength.data(), potential_data.cell_cell_adhesion_strength.data(),
			potential_data.relative_maximum_adhesion_distance.data(), data.agent_type_indices.data(),
			potential_data.cell_adhesion_affinities.data(), data.is_movable.data(), data.neighbors.data());
}

void update_spring_attachments_internal(
	index_t agents_count, real_t time_step, index_t cell_defs_count, const real_t* __restrict__ detachment_rate,
	const real_t* __restrict__ attachment_rate, const real_t* __restrict__ cell_adhesion_affinities,
	const index_t* __restrict__ maximum_number_of_attachments, const index_t* __restrict__ cell_definition_index,
	const std::vector<index_t>* __restrict__ neighbors, std::vector<index_t>* __restrict__ springs)
{
	constexpr index_t erased_spring = -1;

// mark springs for detachment
#pragma omp for
	for (index_t this_cell_index = 0; this_cell_index < agents_count; this_cell_index++)
	{
		for (index_t j = 0; j < (index_t)springs[this_cell_index].size(); j++)
		{
			if (random::instance().uniform() <= detachment_rate[this_cell_index] * time_step)
			{
#pragma omp critical
				{
					const index_t other_cell_index = springs[this_cell_index][j];

					if (other_cell_index != erased_spring)
					{
						springs[this_cell_index][j] = erased_spring;

						*std::find(springs[other_cell_index].begin(), springs[other_cell_index].end(),
								   this_cell_index) = erased_spring;
					}
				}
			}
		}
	}

// remove marked springs
#pragma omp for
	for (index_t this_cell_index = 0; this_cell_index < agents_count; this_cell_index++)
	{
		auto it = std::remove(springs[this_cell_index].begin(), springs[this_cell_index].end(), erased_spring);

		springs[this_cell_index].erase(it, springs[this_cell_index].end());
	}

	// attach cells to springs

#pragma omp for
	for (index_t this_cell_index = 0; this_cell_index < agents_count; this_cell_index++)
	{
		for (std::size_t j = 0; j < neighbors[this_cell_index].size(); j++)
		{
			const index_t other_cell_index = neighbors[this_cell_index][j];

			if (other_cell_index < this_cell_index)
				continue;

			const real_t affinity_l =
				cell_adhesion_affinities[this_cell_index * cell_defs_count + cell_definition_index[other_cell_index]];

			const real_t attachment_prob_l = attachment_rate[this_cell_index] * time_step * affinity_l;

			const real_t affinity_r =
				cell_adhesion_affinities[other_cell_index * cell_defs_count + cell_definition_index[this_cell_index]];

			const real_t attachment_prob_r = attachment_rate[other_cell_index] * time_step * affinity_r;

			if (random::instance().uniform() <= attachment_prob_l || random::instance().uniform() <= attachment_prob_r)
			{
#pragma omp critical
				{
					if ((index_t)springs[this_cell_index].size() < maximum_number_of_attachments[this_cell_index]
						&& (index_t)springs[other_cell_index].size() < maximum_number_of_attachments[other_cell_index])
					{
						springs[this_cell_index].push_back(other_cell_index);
						springs[other_cell_index].push_back(this_cell_index);
					}
				}
			}
		}
	}
}


void base_potential_model::attach_detach_springs(mech_environment& me)
{
	auto& data = me.agent_data;
	auto& potential_data = static_cast<base_potential_data&>(*data.potential_data.get());

	update_spring_attachments_internal(
		data.agents_count(), me.timestep, me.agent_types_count, potential_data.detachment_rate.data(),
		potential_data.attachment_rate.data(), potential_data.cell_adhesion_affinities.data(),
		potential_data.maximum_number_of_attachments.data(), data.agent_type_indices.data(), data.neighbors.data(),
		potential_data.springs.data());
}

template <index_t dims>
void spring_contract_function(index_t agents_count, index_t cell_defs_count, real_t* __restrict__ velocity,
							  const index_t* __restrict__ cell_definition_index,
							  const real_t* __restrict__ attachment_elastic_constant,
							  const real_t* __restrict__ cell_adhesion_affinity, const real_t* __restrict__ position,
							  const std::uint8_t* __restrict__ is_movable, std::vector<index_t>* __restrict__ springs)
{
#pragma omp for
	for (index_t this_cell_index = 0; this_cell_index < agents_count; this_cell_index++)
	{
		if (is_movable[this_cell_index] == 0)
			continue;

		for (std::size_t j = 0; j < springs[this_cell_index].size(); j++)
		{
			const index_t other_cell_index = springs[this_cell_index][j];

			const index_t this_cell_def_index = cell_definition_index[this_cell_index];
			const index_t other_cell_def_index = cell_definition_index[other_cell_index];

			const real_t adhesion =
				sqrt(attachment_elastic_constant[this_cell_index] * attachment_elastic_constant[other_cell_index]
					 * cell_adhesion_affinity[this_cell_index * cell_defs_count + other_cell_def_index]
					 * cell_adhesion_affinity[other_cell_index * cell_defs_count + this_cell_def_index]);

			real_t difference[dims];

			potentials_helper<dims>::subtract(difference, position + other_cell_index * dims,
											  position + this_cell_index * dims);

			potentials_helper<dims>::update_velocity(velocity + this_cell_index * dims, difference, adhesion);
		}
	}
}

void base_potential_model::compute_springs_potentials(mech_environment& me)
{
	auto& data = me.agent_data;
	auto& potential_data = static_cast<base_potential_data&>(*data.potential_data.get());

	if (me.m.mesh.dims == 1)
		spring_contract_function<1>(
			data.agents_count(), me.agent_types_count, data.velocity.data(), data.agent_type_indices.data(),
			potential_data.attachment_elastic_constant.data(), potential_data.cell_adhesion_affinities.data(),
			data.bio_agent_data.positions.data(), data.is_movable.data(), potential_data.springs.data());
	else if (me.m.mesh.dims == 2)
		spring_contract_function<2>(
			data.agents_count(), me.agent_types_count, data.velocity.data(), data.agent_type_indices.data(),
			potential_data.attachment_elastic_constant.data(), potential_data.cell_adhesion_affinities.data(),
			data.bio_agent_data.positions.data(), data.is_movable.data(), potential_data.springs.data());
	else if (me.m.mesh.dims == 3)
		spring_contract_function<3>(
			data.agents_count(), me.agent_types_count, data.velocity.data(), data.agent_type_indices.data(),
			potential_data.attachment_elastic_constant.data(), potential_data.cell_adhesion_affinities.data(),
			data.bio_agent_data.positions.data(), data.is_movable.data(), potential_data.springs.data());
}

void base_potential_model::update_velocities(mech_environment& me)
{
	compute_agents_potentials(me);
	attach_detach_springs(me);
	compute_springs_potentials(me);
}

template <index_t dims>
void update_positions_internal(index_t agents_count, real_t time_step, real_t* __restrict__ position,
							   real_t* __restrict__ velocity, real_t* __restrict__ previous_velocity,
							   const std::uint8_t* __restrict__ is_movable)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		if (!is_movable[i])
			continue;

		const real_t factor = time_step * 1.5;
		const real_t previous_factor = time_step * -0.5;

		for (index_t d = 0; d < dims; d++)
		{
			position[i * dims + d] +=
				velocity[i * dims + d] * factor + previous_velocity[i * dims + d] * previous_factor;

			previous_velocity[i * dims + d] = velocity[i * dims + d];
			velocity[i * dims + d] = 0;
		}
	}
}

void base_potential_model::update_positions(mech_environment& me)
{
	auto& data = me.agent_data;
	auto& potential_data = static_cast<base_potential_data&>(*data.potential_data.get());

	if (me.m.mesh.dims == 1)
		update_positions_internal<1>(data.agents_count(), me.timestep, data.bio_agent_data.positions.data(),
									 data.velocity.data(), potential_data.previous_velocity.data(),
									 data.is_movable.data());
	else if (me.m.mesh.dims == 2)
		update_positions_internal<2>(data.agents_count(), me.timestep, data.bio_agent_data.positions.data(),
									 data.velocity.data(), potential_data.previous_velocity.data(),
									 data.is_movable.data());
	else if (me.m.mesh.dims == 3)
		update_positions_internal<3>(data.agents_count(), me.timestep, data.bio_agent_data.positions.data(),
									 data.velocity.data(), potential_data.previous_velocity.data(),
									 data.is_movable.data());
}
