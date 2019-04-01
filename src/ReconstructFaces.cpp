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
			else
			{
				block_array_4_comparaison_[i][j] = 0;
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
	if(block_array_4_comparaison_[first_block_id][second_block_id]==0)
	{
		std::vector<int>::iterator it;

		std::vector<int> common_nodes_vector(n_nodes_in_first_block+n_nodes_in_second_block);

		std::sort(first_block_node_array,first_block_node_array+n_nodes_in_first_block);
		std::sort(second_block_node_array,second_block_node_array+n_nodes_in_second_block);

		it = std::set_intersection( first_block_node_array, first_block_node_array + n_nodes_in_first_block, second_block_node_array, second_block_node_array + n_nodes_in_second_block, common_nodes_vector.begin() );
		common_nodes_vector.resize(it-common_nodes_vector.begin());

		if(common_nodes_vector.size()!=0)
		{
			std::vector<int> block_ids={first_block_id,second_block_id};
			connexionVector_.push_back(block_ids);

			commonNodesVector_.push_back(common_nodes_vector);


			block_array_4_comparaison_[first_block_id][second_block_id] =1;
			block_array_4_comparaison_[second_block_id][first_block_id] =1;

		}
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
	std::vector<std::vector<int>> all_face_2_nodes_connectivity;

	int first_block_id;
	int second_block_id;

	commonFacesVector_ = new std::vector<std::vector<int>>[n_connexions];

	for(int conn = 0; conn < n_connexions; conn++)
	{
		std::vector<std::vector<int>> common_faces_vector_temp;

		first_block_id = connexionVector_[conn][0];
		second_block_id = connexionVector_[conn][1];

		common_nodes_vector = commonNodesVector_[conn];
		n_common_nodes = common_nodes_vector.size();

		for(int i = 0; i<n_common_nodes;i++)
		{
			node_id = common_nodes_vector[i];
			node_2_cells_connectivity = global_nodes_vector[0][node_id];
			n_cells_in_node = node_2_cells_connectivity.size();

			for(int j = 0; j<n_cells_in_node;j++)
			{
				cell_id = node_2_cells_connectivity[j];

				if(cell_flag_4_face_reconstruction_[cell_id]==0)
				{
					cell_2_nodes_connectivity = global_cells_vector[0][cell_id];
					n_nodes_in_cell = cell_2_nodes_connectivity.size();
					all_face_2_nodes_connectivity = CreateFaces(cell_2_nodes_connectivity);

					for(int k = 0; k<all_face_2_nodes_connectivity.size();k++)
					{
						std::vector<int> face_2_nodes_connectivity = all_face_2_nodes_connectivity [k];

						std::vector<int>::iterator it;
						std::vector<int> face_in_boundary(n_nodes_in_cell+n_common_nodes);
						std::sort(face_2_nodes_connectivity.begin(), face_2_nodes_connectivity.end());
						std::sort(common_nodes_vector.begin(),common_nodes_vector.end());

						it = std::set_intersection( common_nodes_vector.begin(), common_nodes_vector.end(), face_2_nodes_connectivity.begin(), face_2_nodes_connectivity.end(), face_in_boundary.begin() );
						face_in_boundary.resize(it - face_in_boundary.begin());

						// for(int l=0;l<face_in_boundary.size();l++)
						// {
						// 	std::cout << face_in_boundary[l] << '\t';
						// }
						// std::cout << '\n';

						if(face_in_boundary.size()>=3)
						{
							common_faces_vector_temp.push_back(all_face_2_nodes_connectivity [k]);
						}
					}
				}
				cell_flag_4_face_reconstruction_[cell_id] = 1;
			}
		}


		for(int i=0; i<common_faces_vector_temp.size();i++)
		{
			std::vector<int> face_checked = common_faces_vector_temp[i];

			if(face_checked.size()!=0)
			{
				std::sort(face_checked.begin(),face_checked.end());

				for(int j=0; j<common_faces_vector_temp.size();j++)
				{
					if(i!=j)
					{
						std::vector<int> face_in_common_face_vector_temp= common_faces_vector_temp[j];

						if(face_in_common_face_vector_temp.size()!=0)
						{
							std::sort(face_in_common_face_vector_temp.begin(),face_in_common_face_vector_temp.end());

							int same_node_count=0;
							for(int node = 0; node<3;node++)
							{
								if(face_checked[node]==face_in_common_face_vector_temp[node])
								{
									same_node_count+=1;
								}
							}
							if(same_node_count==3)
							{
								commonFacesVector_[conn].push_back(common_faces_vector_temp[i]);
								common_faces_vector_temp[i].clear();
							}
						}
					}
				}
			}
		}
	}
}

std::vector<std::vector<int>> ReconstructFaces::CreateFaces(std::vector<int> cell2NodesConnectivity)
{
	int numberOfNodes = cell2NodesConnectivity.size();
	std::vector<std::vector<int>> face_2_nodes_connectivity_local;
	std::vector<std::vector<int>> face_2_nodes_connectivity;

	// Tetrahedral
	if (numberOfNodes == 4)
	{
		face_2_nodes_connectivity_local = {{0,2,1},{0,1,3},{1,2,3},{2,0,3}};
	}
	// Hexahedral
	else if (numberOfNodes == 8)
	{
		face_2_nodes_connectivity_local = {{1,0,3,2},{1,2,6,5},{2,3,7,6},{3,0,4,7},{0,1,5,4},{4,5,6,7}};
	}
	// Pyramid
	else if (numberOfNodes == 5)
	{
		face_2_nodes_connectivity_local = {{0,3,2,1},{0,1,4},{1,2,4},{2,3,4},{3,0,4}};
	}

	for(int i=0;i<face_2_nodes_connectivity_local.size();i++)
	{
		std::vector<int> face_2_nodes_connectivity_temp;
		for(int j=0;j<face_2_nodes_connectivity_local[i].size();j++)
		{
			int index_2_push_back = face_2_nodes_connectivity_local[i][j];

			face_2_nodes_connectivity_temp.push_back(cell2NodesConnectivity[index_2_push_back]);
		}
		face_2_nodes_connectivity.push_back(face_2_nodes_connectivity_temp);
	}

	return face_2_nodes_connectivity;
}
