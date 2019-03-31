#ifndef MESH_PARTITIONING_HEAD_RECONSTRUCTFACES_H
#define MESH_PARTITIONING_HEAD_RECONSTRUCTFACES_H

#include <iostream>     // std::cout
#include <algorithm>    // std::set_intersection, std::sort
#include <vector>

class ReconstructFaces
{
public:
	int n_blocks_;
	int global_n_elements_;
	int** block_array_4_comparaison_;
	int* cell_flag_4_face_reconstruction_;
	std::vector<std::vector<int>> commonNodesVector_;
	std::vector<std::vector<int>> connexionVector_;
	std::vector<std::vector<int>> commonFacesVector_;


	ReconstructFaces(int n_blocks, int global_n_elements);
	~ReconstructFaces();

	void CompareArraysOfGlobalNodes( int first_block_id, int* first_block_node_array, int n_nodes_in_first_block, int second_block_id, int* second_block_node_array, int n_nodes_in_second_block);
	void FindElementsInConnexion(std::vector<int>** global_cells_vector, std::vector<int>** global_nodes_vector);
	std::vector<std::vector<int>> CreateFaces(std::vector<int> cell2NodesConnectivity);

};

#endif
