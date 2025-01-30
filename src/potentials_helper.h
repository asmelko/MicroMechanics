#pragma once

#include <cmath>

#include <BioFVM/types.h>

namespace micromech {

template <biofvm::index_t dims>
struct potentials_helper
{};

template <>
struct potentials_helper<1>
{
	static constexpr biofvm::real_t distance(const biofvm::real_t* __restrict__ lhs,
											 const biofvm::real_t* __restrict__ rhs)
	{
		return std::abs(lhs[0] - rhs[0]);
	}

	static constexpr biofvm::real_t difference_and_distance(const biofvm::real_t* __restrict__ lhs,
															const biofvm::real_t* __restrict__ rhs,
															biofvm::real_t* __restrict__ difference)
	{
		difference[0] = lhs[0] - rhs[0];

		return std::abs(difference[0]);
	}

	static constexpr void update_velocity(biofvm::real_t* __restrict__ velocity,
										  const biofvm::real_t* __restrict__ difference, const biofvm::real_t force)
	{
		velocity[0] += force * difference[0];
	}

	static constexpr void subtract(biofvm::real_t* __restrict__ dst, const biofvm::real_t* __restrict__ lhs,
								   const biofvm::real_t* __restrict__ rhs)
	{
		dst[0] = lhs[0] - rhs[0];
	}

	static constexpr void add(biofvm::real_t* __restrict__ lhs, const biofvm::real_t* __restrict__ rhs)
	{
		lhs[0] += rhs[0];
	}
};

template <>
struct potentials_helper<2>
{
	static constexpr biofvm::real_t distance(const biofvm::real_t* __restrict__ lhs,
											 const biofvm::real_t* __restrict__ rhs)
	{
		return std::sqrt((lhs[0] - rhs[0]) * (lhs[0] - rhs[0]) + (lhs[1] - rhs[1]) * (lhs[1] - rhs[1]));
	}

	static constexpr biofvm::real_t difference_and_distance(const biofvm::real_t* __restrict__ lhs,
															const biofvm::real_t* __restrict__ rhs,
															biofvm::real_t* __restrict__ difference)
	{
		difference[0] = lhs[0] - rhs[0];
		difference[1] = lhs[1] - rhs[1];

		return std::sqrt(difference[0] * difference[0] + difference[1] * difference[1]);
	}

	static constexpr void update_velocity(biofvm::real_t* __restrict__ velocity,
										  const biofvm::real_t* __restrict__ difference, const biofvm::real_t force)
	{
		velocity[0] += force * difference[0];
		velocity[1] += force * difference[1];
	}

	static constexpr void subtract(biofvm::real_t* __restrict__ dst, const biofvm::real_t* __restrict__ lhs,
								   const biofvm::real_t* __restrict__ rhs)
	{
		dst[0] = lhs[0] - rhs[0];
		dst[1] = lhs[1] - rhs[1];
	}

	static constexpr void add(biofvm::real_t* __restrict__ lhs, const biofvm::real_t* __restrict__ rhs)
	{
		lhs[0] += rhs[0];
		lhs[1] += rhs[1];
	}
};

template <>
struct potentials_helper<3>
{
	static constexpr biofvm::real_t distance(const biofvm::real_t* __restrict__ lhs,
											 const biofvm::real_t* __restrict__ rhs)
	{
		return std::sqrt((lhs[0] - rhs[0]) * (lhs[0] - rhs[0]) + (lhs[1] - rhs[1]) * (lhs[1] - rhs[1])
						 + (lhs[2] - rhs[2]) * (lhs[2] - rhs[2]));
	}

	static constexpr biofvm::real_t difference_and_distance(const biofvm::real_t* __restrict__ lhs,
															const biofvm::real_t* __restrict__ rhs,
															biofvm::real_t* __restrict__ difference)
	{
		difference[0] = lhs[0] - rhs[0];
		difference[1] = lhs[1] - rhs[1];
		difference[2] = lhs[2] - rhs[2];

		return std::sqrt(difference[0] * difference[0] + difference[1] * difference[1] + difference[2] * difference[2]);
	}

	static constexpr void update_velocity(biofvm::real_t* __restrict__ velocity,
										  const biofvm::real_t* __restrict__ difference, const biofvm::real_t force)
	{
		velocity[0] += force * difference[0];
		velocity[1] += force * difference[1];
		velocity[2] += force * difference[2];
	}

	static constexpr void subtract(biofvm::real_t* __restrict__ dst, const biofvm::real_t* __restrict__ lhs,
								   const biofvm::real_t* __restrict__ rhs)
	{
		dst[0] = lhs[0] - rhs[0];
		dst[1] = lhs[1] - rhs[1];
		dst[2] = lhs[2] - rhs[2];
	}

	static constexpr void add(biofvm::real_t* __restrict__ lhs, const biofvm::real_t* __restrict__ rhs)
	{
		lhs[0] += rhs[0];
		lhs[1] += rhs[1];
		lhs[2] += rhs[2];
	}
};

} // namespace micromech
