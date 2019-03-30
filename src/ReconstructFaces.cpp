#include "ReconstructFaces.h"

ReconstructFaces::ReconstructFaces(int n_blocks, int global_n_elements):n_blocks_(n_blocks), global_n_elements_(global_n_elements), block_array_4_comparaison_(nullptr), cell_flag_4_face_reconstruction_(nullptr)
{
	block_array_4_comparaison_ = new int*[n_blocks_]();

	cell_flag_4_face_reconstruction_ = new int[global_n_elements_]();

	for(int i = 0; i< n_blocks_; i++)
	{
		block_array_4_comparaison_[i] = new int[n_blocks_]();
		for(int j = 0; j < n_blocks_; j++)
		{
			if(i==j)
			{
				block_array_4_comparaison_[i][j] = 1;
			}
		}
	}

}

ReconstructFaces::~ReconstructFaces()
{
	if(block_array_4_comparaison_!=nullptr)
	{
		for(int i = 0; i< n_blocks_; i++)
		{
			if(block_array_4_comparaison_[i]!=nullptr)
			{
				delete[] block_array_4_comparaison_[i];
				block_array_4_comparaison_[i] = nullptr;
			}
		}
		delete[] block_array_4_comparaison_;
		block_array_4_comparaison_=nullptr;
	}

	if(cell_flag_4_face_reconstruction_!=nullptr)
	{
		delete[] cell_flag_4_face_reconstruction_;
		cell_flag_4_face_reconstruction_ = nullptr;
	}

}



void ReconstructFaces::CompareArraysOfGlobalNodes( int first_block_id, int* first_block_node_array, int n_nodes_in_first_block, int second_block_id, int* second_block_node_array, int n_nodes_in_second_block)
{

	std::vector<int>::iterator it;

	std::vector<int> common_nodes_vector(n_nodes_in_first_block+n_nodes_in_second_block);

	std::sort(first_block_node_array,first_block_node_array+n_nodes_in_first_block);
	std::sort(second_block_node_array,second_block_node_array+n_nodes_in_second_block);

	it = std::set_intersection( first_block_node_array, first_block_node_array + n_nodes_in_first_block, second_block_node_array, second_block_node_array + n_nodes_in_second_block, common_nodes_vector.begin() );

	if(!common_nodes_vector.empty())
	{
		common_nodes_vector.resize(it-common_nodes_vector.begin());
		std::vector<int> block_ids={first_block_id,second_block_id};
		connexionVector_.push_back(block_ids);

		commonNodesVector_.push_back(common_nodes_vector);


		block_array_4_comparaison_[first_block_id][second_block_id] =1;
		block_array_4_comparaison_[second_block_id][first_block_id] =1;

	}
}

void ReconstructFaces::FindElementsInConnexion(std::vector<int>** global_cells_vector, std::vector<int>** global_nodes_vector)
{
	int node_id;
	int cell_id;
	int n_cells_in_node;
	int n_nodes_in_cell;
	int n_common_nodes;

	int n_connexions = connexionVector_.size();
	std::vector<int> common_nodes_vector;
	std::vector<int> node_2_cells_connectivity;
	std::vector<int> cell_2_nodes_connectivity;



	int block_id;

	for(int conn = 0; conn < n_connexions; conn++)
	{
		common_nodes_vector = commonNodesVector_[conn];
		n_common_nodes = common_nodes_vector.size();

		for(int blockI; blockI<2;blockI++)
		{
			block_id = connexionVector_[conn][blockI];
			for(int i=0; i<n_common_nodes; i++)
			{
				node_id = common_nodes_vector[i];
				node_2_cells_connectivity = global_nodes_vector[0][node_id];

				n_cells_in_node = node_2_cells_connectivity.size();

				for(int j=0; j<n_cells_in_node; j++)
				{
					cell_id = node_2_cells_connectivity[j];

					cell_2_nodes_connectivity = global_cells_vector[0][cell_id];

					n_nodes_in_cell = cell_2_nodes_connectivity.size();

					if(cell_flag_4_face_reconstruction_[cell_id]==0)
					{
						std::vector<int>::iterator it;
						std::vector<int> face_2_nodes_connectivity(4);

						std::sort(cell_2_nodes_connectivity.begin(), cell_2_nodes_connectivity.end());

						it = std::set_intersection( common_nodes_vector.begin(), common_nodes_vector.end(), cell_2_nodes_connectivity.begin(), cell_2_nodes_connectivity.begin(), face_2_nodes_connectivity.begin() );

						if(!face_2_nodes_connectivity.empty())
						{
							face_2_nodes_connectivity.resize(it - face_2_nodes_connectivity.begin());
							commonFacesVector_.push_back(face_2_nodes_connectivity);
							cell_flag_4_face_reconstruction_[cell_id] = 1;


						}
					}
				}
			}

		}

	}






}
