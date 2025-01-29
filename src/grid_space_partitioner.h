#pragma once

#include <atomic>
#include <memory>
#include <vector>

#include <BioFVM/mesh.h>

#include "BioFVM/types.h"

namespace micromech {

class grid_space_partitioner
{
	biofvm::cartesian_mesh partitioning_mesh_;

	std::unique_ptr<std::atomic<biofvm::index_t>[]> agents_in_voxels_sizes_;
	std::unique_ptr<std::vector<biofvm::index_t>[]> agents_in_voxels_;

	template <biofvm::index_t dims>
	biofvm::index_t get_mesh_index(biofvm::point_t<biofvm::index_t, 3> point) const;

	biofvm::index_t get_mesh_index(const biofvm::real_t* position) const;

public:
	grid_space_partitioner(biofvm::index_t voxel_size, const biofvm::cartesian_mesh& microenv_mesh);

	void update_partitioning(const biofvm::real_t* positions, biofvm::index_t agents_count);

	template <biofvm::index_t dims, typename func_t>
	void for_each_in_neighborhood(const biofvm::real_t* agent_position, biofvm::index_t i, func_t f)
	{
		auto position = partitioning_mesh_.voxel_position<dims>(agent_position);
		for (biofvm::index_t z = -1; z <= 1; z++)
		{
			if (position[2] + z >= partitioning_mesh_.grid_shape[2] || position[2] + z < 0)
				continue;

			for (biofvm::index_t y = -1; y <= 1; y++)
			{
				if (position[1] + y >= partitioning_mesh_.grid_shape[1] || position[1] + y < 0)
					continue;

				for (biofvm::index_t x = -1; x <= 1; x++)
				{
					if (position[0] + x >= partitioning_mesh_.grid_shape[0] || position[0] + x < 0)
						continue;

					biofvm::index_t voxel_index;

					if constexpr (dims == 1)
					{
						voxel_index = get_mesh_index<1>({ position[0] + x, 0, 0 });
					}
					else if constexpr (dims == 2)
					{
						voxel_index = get_mesh_index<2>({ position[0] + x, position[1] + y, 0 });
					}
					else
					{
						voxel_index = get_mesh_index<3>({ position[0] + x, position[1] + y, position[2] + z });
					}

					for (auto& cell_idx : agents_in_voxels_[voxel_index])
					{
						if (i != cell_idx)
							f(cell_idx);
					}
				}
			}
		}
	}
};

} // namespace micromech
