#ifndef MESH_PARTITIONING_HEAD_RECONSTRUCTFACES_H
#define MESH_PARTITIONING_HEAD_RECONSTRUCTFACES_H

#include <iostream>     // std::cout
#include <algorithm>    // std::set_intersection, std::sort
#include <vector>
#include "MetisBoundary.h"


class ReconstructFaces
{
public:
	int n_blocks_;
	int global_n_elements_;
	int** block_array_4_comparaison_;
	int* cell_flag_4_face_reconstruction_;
	std::vector<std::vector<int>> commonNodesVector_;
	std::vector<std::vector<int>> connexionVector_;
	std::vector<std::vector<int>>* commonFacesVector_;
	std::vector<std::vector<int>>* commonCellsVector_;
	std::vector<std::vector<int>>* boundaryFaces_;
	std::vector<int>* boundaryElements_;


	ReconstructFaces(int n_blocks, int global_n_elements);
	~ReconstructFaces();

	void CompareArraysOfGlobalNodes( int first_block_id, int* first_block_node_array, int n_nodes_in_first_block, int second_block_id, int* second_block_node_array, int n_nodes_in_second_block);
	void FindElementsInConnexion(std::vector<int>** global_cells_vector, std::vector<int>** global_nodes_vector);
	std::vector<std::vector<int>> CreateFaces(std::vector<int> cell2NodesConnectivity);
	void FindBoundaryFaces(int* nElements,std::vector<int>** cell_2_nodes_connectivity,std::vector<int>** node_2_cells_connectivity);
	void FindConnexionFaces(int** local_2_global_nodes,std::vector<int>* global_2_local_elements,std::vector<int>** global_node_2_cells_connectivity);
	void RemovePhysicalBoundaries(int** local_2_global_nodes, MetisBoundary* single_block_metisboundary);

};

#endif
