#include "base_motility_model.h"

#include "base_motility_data.h"
#include "potentials_helper.h"
#include "random.h"

using namespace biofvm;
using namespace micromech;

constexpr biofvm::real_t zero_threshold = 1e-16;

template <index_t dims>
struct motility_helper
{};

template <>
struct motility_helper<1>
{
	static void random_walk(bool, biofvm::real_t* __restrict__ walk)
	{
		biofvm::real_t rand = random::instance().uniform();
		walk[0] = rand < 0.5 ? -1 : 1;
	}

	static constexpr void update_motility_vector(biofvm::real_t* __restrict__ motility_vector,
												 const biofvm::real_t* __restrict__ walk,
												 const biofvm::real_t* __restrict__ migration_bias_direction,
												 const biofvm::real_t migration_bias)
	{
		motility_vector[0] = (1 - migration_bias) * walk[0] + migration_bias * migration_bias_direction[0];
	}

	static constexpr void normalize_and_scale(biofvm::real_t* __restrict__ vector, biofvm::real_t scale)
	{
		biofvm::real_t length = std::abs(vector[0]);

		vector[0] = length > zero_threshold ? vector[0] * scale / length : 0;
	}
};

template <>
struct motility_helper<2>
{
	static void random_walk(bool, biofvm::real_t* __restrict__ walk)
	{
		biofvm::real_t theta = random::instance().uniform(0, 2 * std::numbers::pi_v<biofvm::real_t>);
		walk[0] = std::cos(theta);
		walk[1] = std::sin(theta);
	}

	static constexpr void update_motility_vector(biofvm::real_t* __restrict__ motility_vector,
												 const biofvm::real_t* __restrict__ walk,
												 const biofvm::real_t* __restrict__ migration_bias_direction,
												 const biofvm::real_t migration_bias)
	{
		motility_vector[0] = (1 - migration_bias) * walk[0] + migration_bias * migration_bias_direction[0];
		motility_vector[1] = (1 - migration_bias) * walk[1] + migration_bias * migration_bias_direction[1];
	}

	static constexpr void normalize_and_scale(biofvm::real_t* __restrict__ vector, biofvm::real_t scale)
	{
		biofvm::real_t length = std::sqrt(vector[0] * vector[0] + vector[1] * vector[1]);

		vector[0] = length > zero_threshold ? vector[0] * scale / length : 0;
		vector[1] = length > zero_threshold ? vector[1] * scale / length : 0;
	}
};

template <>
struct motility_helper<3>
{
	static void random_walk(bool restrict_to_2d, biofvm::real_t* __restrict__ walk)
	{
		if (restrict_to_2d)
		{
			motility_helper<2>::random_walk(true, walk);
			walk[2] = 0;
		}
		else
		{
			const biofvm::real_t theta = random::instance().uniform(0, 2 * std::numbers::pi_v<biofvm::real_t>);
			const biofvm::real_t z = random::instance().uniform(-1, 1);
			const biofvm::real_t r = std::sqrt(1 - z * z);

			walk[0] = std::cos(theta) * r;
			walk[1] = std::sin(theta) * r;
			walk[2] = z;
		}
	}

	static constexpr void update_motility_vector(biofvm::real_t* __restrict__ motility_vector,
												 const biofvm::real_t* __restrict__ walk,
												 const biofvm::real_t* __restrict__ migration_bias_direction,
												 const biofvm::real_t migration_bias)
	{
		motility_vector[0] = (1 - migration_bias) * walk[0] + migration_bias * migration_bias_direction[0];
		motility_vector[1] = (1 - migration_bias) * walk[1] + migration_bias * migration_bias_direction[1];
		motility_vector[2] = (1 - migration_bias) * walk[2] + migration_bias * migration_bias_direction[2];
	}

	static constexpr void normalize_and_scale(biofvm::real_t* __restrict__ vector, biofvm::real_t scale)
	{
		biofvm::real_t length = std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);

		vector[0] = length > zero_threshold ? vector[0] * scale / length : 0;
		vector[1] = length > zero_threshold ? vector[1] * scale / length : 0;
		vector[2] = length > zero_threshold ? vector[2] * scale / length : 0;
	}
};

base_motility_model::base_motility_model(mech_environment& me)
{
	if (dynamic_cast<base_motility_data*>(me.agent_data.motility_data.get()) == nullptr)
	{
		me.agent_data.motility_data = std::make_unique<base_motility_data>(me);
	}
}

template <index_t dims>
void update_motility_internal(
	index_t agents_count, real_t time_step, real_t* __restrict__ motility_vector, real_t* __restrict__ velocity,
	const real_t* __restrict__ persistence_time, const real_t* __restrict__ migration_bias,
	real_t* __restrict__ migration_bias_direction, const std::uint8_t* __restrict__ restrict_to_2d,
	const std::uint8_t* __restrict__ is_motile, const real_t* __restrict__ migration_speed,
	const base_motility_data::direction_update_func* __restrict__ update_migration_bias_direction_f)
{
#pragma omp for
	for (index_t i = 0; i < agents_count; i++)
	{
		if (is_motile[i] == 0)
			continue;

		if (random::instance().uniform() < time_step / persistence_time[i])
		{
			real_t random_walk[dims];

			motility_helper<dims>::random_walk(restrict_to_2d, random_walk);

			if (update_migration_bias_direction_f[i] != nullptr)
			{
				update_migration_bias_direction_f[i](migration_bias_direction + i * dims);
			}

			motility_helper<dims>::update_motility_vector(motility_vector + i * dims, random_walk,
														  migration_bias_direction + i * dims, migration_bias[i]);

			motility_helper<dims>::normalize_and_scale(motility_vector + i * dims, migration_speed[i]);
		}

		potentials_helper<dims>::add(velocity + i * dims, motility_vector + i * dims);
	}
}

void base_motility_model::update_motility_velocities(mech_environment& me)
{
	auto& data = me.agent_data;
	auto& motility_data = static_cast<base_motility_data&>(*data.motility_data.get());


	if (me.m.mesh.dims == 1)
		update_motility_internal<1>(data.agents_count(), me.timestep, motility_data.motility_vector.data(),
									data.velocity.data(), motility_data.persistence_time.data(),
									motility_data.migration_bias.data(), motility_data.migration_bias_direction.data(),
									motility_data.restrict_to_2d.data(), motility_data.is_motile.data(),
									motility_data.migration_speed.data(),
									motility_data.update_migration_bias_direction.data());
	else if (me.m.mesh.dims == 2)
		update_motility_internal<2>(data.agents_count(), me.timestep, motility_data.motility_vector.data(),
									data.velocity.data(), motility_data.persistence_time.data(),
									motility_data.migration_bias.data(), motility_data.migration_bias_direction.data(),
									motility_data.restrict_to_2d.data(), motility_data.is_motile.data(),
									motility_data.migration_speed.data(),
									motility_data.update_migration_bias_direction.data());
	else if (me.m.mesh.dims == 3)
		update_motility_internal<3>(data.agents_count(), me.timestep, motility_data.motility_vector.data(),
									data.velocity.data(), motility_data.persistence_time.data(),
									motility_data.migration_bias.data(), motility_data.migration_bias_direction.data(),
									motility_data.restrict_to_2d.data(), motility_data.is_motile.data(),
									motility_data.migration_speed.data(),
									motility_data.update_migration_bias_direction.data());
}
