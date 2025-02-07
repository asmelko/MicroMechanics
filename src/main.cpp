#include <chrono>
#include <iostream>
#include <random>

#include <BioFVM/microenvironment.h>

#include "BioFVM/types.h"
#include "base_membrane_data.h"
#include "base_motility_data.h"
#include "base_motility_model.h"
#include "base_potential_data.h"
#include "base_potential_model.h"
#include "base_vtk_serializer.h"
#include "base_wall_membrane_model.h"
#include "grid_space_partitioner.h"
#include "mech_environment.h"

using namespace biofvm;
using namespace micromech;

using data_setup_func_t = std::function<void(index_t, mech_environment&)>;

void setup_base_membrane_data(index_t i, mech_environment& me)
{
	auto membrane_data = dynamic_cast<base_membrane_data*>(me.agent_data.membrane_data.get());

	membrane_data->cell_BM_repulsion_strength[i] = 1;
}

void setup_base_motility_data(index_t i, mech_environment& me)
{
	auto motility_data = dynamic_cast<base_motility_data*>(me.agent_data.motility_data.get());

	motility_data->is_motile[i] = false;
}

void setup_base_potential_data(index_t i, mech_environment& me)
{
	auto potential_data = dynamic_cast<base_potential_data*>(me.agent_data.potential_data.get());

	potential_data->cell_cell_adhesion_strength[i] = 1;
	potential_data->cell_cell_repulsion_strength[i] = 10;
	potential_data->cell_adhesion_affinities[i] = 3;
	potential_data->relative_maximum_adhesion_distance[i] = 1;
	potential_data->maximum_number_of_attachments[i] = 1;
	potential_data->attachment_elastic_constant[i] = 1;
	potential_data->attachment_rate[i] = 0;
	potential_data->detachment_rate[i] = 0;
}

void make_agents(std::size_t count, mech_environment& me, data_setup_func_t&& setup_membrane_data,
				 data_setup_func_t&& setup_motility_data, data_setup_func_t&& setup_potential_data)
{
	std::uniform_real_distribution<real_t> distr_x(me.m.mesh.bounding_box_mins[0] + 200,
												   me.m.mesh.bounding_box_maxs[0] - 200);
	std::uniform_real_distribution<real_t> distr_y(me.m.mesh.bounding_box_mins[1] + 200,
												   me.m.mesh.bounding_box_maxs[1] - 200);
	std::uniform_real_distribution<real_t> distr_z(me.m.mesh.bounding_box_mins[2] + 200,
												   me.m.mesh.bounding_box_maxs[2] - 200);

	std::mt19937 gen;

	for (std::size_t i = 0; i < count; ++i)
	{
		me.agent_data.add();

		me.agent_data.radius[i] = 10;
		me.agent_data.is_movable[i] = true;
		me.agent_data.agent_type_indices[i] = 0;

		for (index_t dim = 0; dim < me.m.mesh.dims; ++dim)
		{
			me.agent_data.bio_agent_data.positions[i * me.m.mesh.dims + dim] =
				std::vector { distr_x, distr_y, distr_z }[dim](gen);
		}

		setup_membrane_data(i, me);
		setup_motility_data(i, me);
		setup_potential_data(i, me);
	}
}

int main()
{
	cartesian_mesh mesh(2, { 0, 0, 0 }, { 1000, 1000, 0 }, { 20, 20, 20 });

	real_t diffusion_time_step = 1;
	index_t substrates_count = 4;
	auto initial_conds = std::make_unique<real_t[]>(substrates_count);

	microenvironment m(mesh, substrates_count, diffusion_time_step, initial_conds.get());

	real_t mech_time_step = 1;
	index_t agent_types_count = 4;

	mech_environment me(m, mech_time_step, agent_types_count);

	me.membrane_m = std::make_unique<base_wall_membrane_model>(me);

	grid_space_partitioner partitioner(20, mesh);
	me.potential_m = std::make_unique<base_potential_model>(partitioner, me);

	me.motility_m = std::make_unique<base_motility_model>(me);

	size_t agents_count = 2000;
	make_agents(agents_count, me, setup_base_membrane_data, setup_base_motility_data, setup_base_potential_data);

	base_vtk_serializer vtk_serializer("output");

#pragma omp parallel
	for (index_t i = 0; i < 500; i++)
	{
		std::size_t partition_duration, membrane_duration, motility_duration, neighbors_duration, velocities_duration,
			positions_duration;

		{
			auto start = std::chrono::high_resolution_clock::now();

			partitioner.update_partitioning(me.agent_data.bio_agent_data.positions.data(), agents_count);

			auto end = std::chrono::high_resolution_clock::now();

			partition_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}

		{
			auto start = std::chrono::high_resolution_clock::now();

			me.membrane_m->compute_basement_membrane_interactions(me);

			auto end = std::chrono::high_resolution_clock::now();

			membrane_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}

		{
			auto start = std::chrono::high_resolution_clock::now();

			me.motility_m->update_motility_velocities(me);

			auto end = std::chrono::high_resolution_clock::now();

			motility_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}

		{
			auto start = std::chrono::high_resolution_clock::now();

			me.potential_m->update_neighbors(me);

			auto end = std::chrono::high_resolution_clock::now();

			neighbors_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}

		{
			auto start = std::chrono::high_resolution_clock::now();

			me.potential_m->update_velocities(me);

			auto end = std::chrono::high_resolution_clock::now();

			velocities_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}

		{
			auto start = std::chrono::high_resolution_clock::now();

			me.potential_m->update_positions(me);

			auto end = std::chrono::high_resolution_clock::now();

			positions_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}

#pragma omp barrier
#pragma omp master
		{
			std::size_t serialization_duration;
			{
				auto start = std::chrono::high_resolution_clock::now();

				vtk_serializer.serialize_one_timestep(me);

				auto end = std::chrono::high_resolution_clock::now();

				serialization_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			}

			std::cout << "Partition time: " << partition_duration << " ms,\t Membrane time: " << membrane_duration
					  << " ms,\t Motility time: " << motility_duration
					  << " ms,\t Neighbors time: " << neighbors_duration
					  << " ms,\t Velocities time: " << velocities_duration
					  << " ms,\t Positions time: " << positions_duration
					  << " ms,\t Serialization time: " << serialization_duration << " ms" << std::endl;
		}
#pragma omp barrier
	}

	return 0;
}
