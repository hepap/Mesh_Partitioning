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

		int* first_block_node_array_2_sort = new int[n_nodes_in_first_block];

		for(int i=0;i<n_nodes_in_first_block;i++)
		{
			first_block_node_array_2_sort[i]= first_block_node_array[i];
		}

		int* second_block_node_array_2_sort = new int[n_nodes_in_second_block];
		for(int i=0;i<n_nodes_in_second_block;i++)
		{
			second_block_node_array_2_sort[i]= second_block_node_array[i];
		}
		std::sort(first_block_node_array_2_sort,first_block_node_array_2_sort+n_nodes_in_first_block);
		std::sort(second_block_node_array_2_sort,second_block_node_array_2_sort+n_nodes_in_second_block);

		it = std::set_intersection( first_block_node_array_2_sort, first_block_node_array_2_sort + n_nodes_in_first_block, second_block_node_array_2_sort, second_block_node_array_2_sort + n_nodes_in_second_block, common_nodes_vector.begin() );
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
	commonCellsVector_ = new std::vector<std::vector<int>>[n_connexions];



	for(int conn = 0; conn < n_connexions; conn++)
	{
		for(int cell_idx = 0; cell_idx<global_n_elements_;cell_idx++)
		{
			cell_flag_4_face_reconstruction_[cell_id]=0;
		}
		std::vector<std::vector<int>> common_faces_vector_temp;
		std::vector<int> common_cells_vector_temp;

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
							common_cells_vector_temp.push_back(cell_id);
						}
					}
				}
				// cell_flag_4_face_reconstruction_[cell_id] = 1;
			}
		}

		std::vector<int> cell_flag(common_cells_vector_temp.size(),0);
		std::vector<int> face_flag(common_faces_vector_temp.size(),0);


		for(int i=0; i<common_faces_vector_temp.size();i++)
		{

			std::vector<int> face_checked = common_faces_vector_temp[i];
			std::vector<int> face_index;
			int face_flag_count =0;

			std::sort(face_checked.begin(),face_checked.end());

			if(face_flag[i]==0)
			{
				for(int j=0; j<common_faces_vector_temp.size();j++)
				{
					if(i!=j)
					{
							std::vector<int> face_in_common_face_vector_temp= common_faces_vector_temp[j];

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
								face_index.push_back(j);
								face_flag_count+=1;
							}
						}
				}
				if(face_flag_count>0)
				{

					commonFacesVector_[conn].push_back(common_faces_vector_temp[i]);
					face_flag[i] =1;
					int common_cell_2_push_back;
					for(int j=0; j<face_flag_count;j++)
					{
						face_flag[face_index[j]]=1;

						if(common_cells_vector_temp[i]!=common_cells_vector_temp[face_index[j]])
						{
							common_cell_2_push_back=common_cells_vector_temp[face_index[j]];
						}
					}
					std::vector<int> common_cells = {common_cells_vector_temp[i],common_cell_2_push_back};
					commonCellsVector_[conn].push_back(common_cells);

				}
		}




		// for(int i=0; i<common_faces_vector_temp.size();i++)
		// {
		//
		// 	if(face_flag[i]!=1)
		// 	{
		// 		std::vector<int> face_checked = common_faces_vector_temp[i];
		//
		// 		std::sort(face_checked.begin(),face_checked.end());
		//
		// 		for(int j=0; j<common_faces_vector_temp.size();j++)
		// 		{
		// 			if(i!=j)
		// 			{
		//
		// 				if(face_flag[j]!=1)
		// 				{
		// 					std::vector<int> face_in_common_face_vector_temp= common_faces_vector_temp[j];
		//
		// 					std::sort(face_in_common_face_vector_temp.begin(),face_in_common_face_vector_temp.end());
		//
		// 					int same_node_count=0;
		// 					for(int node = 0; node<3;node++)
		// 					{
		// 						if(face_checked[node]==face_in_common_face_vector_temp[node])
		// 						{
		// 							same_node_count+=1;
		// 						}
		// 					}
		// 					if(same_node_count==3 && face_flag[j]!=1 && face_flag[i]!=1)
		// 					{
		// 							commonFacesVector_[conn].push_back(common_faces_vector_temp[i]);
		// 							std::vector<int> common_cells = {common_cells_vector_temp[i],common_cells_vector_temp[j]};
		// 							commonCellsVector_[conn].push_back(common_cells);
		// 							// common_faces_vector_temp[i].clear();
		// 							// common_faces_vector_temp[j].clear();
		// 							// face_flag[i]=1;
		// 							face_flag[j]=1;
		//
		// 							cell_flag[i]=1;
		// 							cell_flag[j]=1;
		// 							break;
		// 					}
		// 				}
		// 			}
		// 		}
		// 	}
		// }
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

	for(size_t i=0;i<face_2_nodes_connectivity_local.size();i++)
	{
		std::vector<int> face_2_nodes_connectivity_temp;
		for(size_t j = 0; j < face_2_nodes_connectivity_local[i].size(); j++)
		{
			int index_2_push_back = face_2_nodes_connectivity_local[i][j];

			face_2_nodes_connectivity_temp.push_back(cell2NodesConnectivity[index_2_push_back]);
		}
		face_2_nodes_connectivity.push_back(face_2_nodes_connectivity_temp);
	}

	return face_2_nodes_connectivity;
}

void ReconstructFaces::FindBoundaryFaces(int* nElements,std::vector<int>** cell_2_nodes_connectivity,std::vector<int>** node_2_cells_connectivity)
{
	std::cout<<"In FindBoundaryFaces"<<std::endl;
	boundaryFaces_ = new std::vector<std::vector<int>>[n_blocks_];
	boundaryElements_ = new std::vector<int>[n_blocks_];

	for(int blockI = 0; blockI<n_blocks_;blockI++)
	{
		std::cout<<"======================================================="<<std::endl;

		for(int cellI = 0; cellI<nElements[blockI];cellI++)
		{
			std::vector<int> cell_2_nodes_connectivity_cellI = cell_2_nodes_connectivity[blockI][cellI];
			std::vector<std::vector<int>> faces = CreateFaces(cell_2_nodes_connectivity_cellI);

			for(int faceI = 0; faceI<int(faces.size()); faceI++)
			{
				int node0 = faces[faceI][0];
				int node1 = faces[faceI][1];
				int node2 = faces[faceI][2];

				std::vector<int> node_2_cells_connectivity_node0 = node_2_cells_connectivity[blockI][node0];
				std::vector<int> node_2_cells_connectivity_node1 = node_2_cells_connectivity[blockI][node1];
				std::vector<int> node_2_cells_connectivity_node2 = node_2_cells_connectivity[blockI][node2];

				std::sort(node_2_cells_connectivity_node0.begin(),node_2_cells_connectivity_node0.end());
				std::sort(node_2_cells_connectivity_node1.begin(),node_2_cells_connectivity_node1.end());
				std::sort(node_2_cells_connectivity_node2.begin(),node_2_cells_connectivity_node2.end());

				std::vector<int>::iterator first_it;
				std::vector<int> first_cells_intersection(node_2_cells_connectivity_node0.size()+node_2_cells_connectivity_node1.size());

				first_it = std::set_intersection( node_2_cells_connectivity_node0.begin(), node_2_cells_connectivity_node0.end(), node_2_cells_connectivity_node0.begin(), node_2_cells_connectivity_node0.end(), first_cells_intersection.begin() );
				first_cells_intersection.resize(first_it-first_cells_intersection.begin());

				if(first_cells_intersection.size()!=0)
				{
					std::vector<int>::iterator second_it;
					std::vector<int> second_cells_intersection(node_2_cells_connectivity_node2.size()+first_cells_intersection.size());
					second_it = std::set_intersection( node_2_cells_connectivity_node2.begin(), node_2_cells_connectivity_node2.end(), first_cells_intersection.begin(), first_cells_intersection.end(), second_cells_intersection.begin() );

					second_cells_intersection.resize(second_it - second_cells_intersection.begin());

					if(second_cells_intersection.size()==1)
					{
						std::cout<<"Frontiere "<<second_cells_intersection[0]<<std::endl;
						boundaryFaces_[blockI].push_back(faces[faceI]);
						boundaryElements_[blockI].push_back(second_cells_intersection[0]);
						std::cout<<faces[faceI][0]<<"\t"<<faces[faceI][1]<<"\t"<<faces[faceI][2]<<"\t"<<faces[faceI][3]<<std::endl;


					}
					else if(second_cells_intersection.size()==2)
					{
						// std::cout<<"Normal "<<second_cells_intersection[0]<<"\t"<<second_cells_intersection[1]<<std::endl;

					}
				}
			}
		}
	}
}

void ReconstructFaces::FindConnexionFaces(int** local_2_global_nodes, std::vector<int>* global_2_local_elements,std::vector<int>** global_node_2_cells_connectivity)
{
	commonFacesVector_ = new std::vector<std::vector<int>>[connexionVector_.size()];
	std::vector<std::vector<int>>* common_faces_vector_temp= new std::vector<std::vector<int>>[connexionVector_.size()];
	commonCellsVector_ = new std::vector<std::vector<int>>[connexionVector_.size()];

	for(int connexionI =0; connexionI < int(connexionVector_.size()); connexionI++)
	{

		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" << '\n';

		int count =0;
		std::cout << connexionI << '\n';
		int first_block_id = connexionVector_[connexionI][0];
		int second_block_id = connexionVector_[connexionI][1];

		std::cout<< first_block_id<<std::endl;
		std::vector<int> common_nodes_vector = commonNodesVector_[connexionI];
		int blockI = first_block_id;
		for(int faceI=0; faceI<boundaryFaces_[blockI].size();faceI++)
		{
			if(boundaryFaces_[blockI][faceI].size()==0)
			{
				continue;
			}

			int cell_local = boundaryElements_[blockI][faceI];

			std::vector<int> face_2_nodes_connectivity;
			for(int nodeI=0;nodeI<boundaryFaces_[blockI][faceI].size();nodeI++)
			{
				int local_node = boundaryFaces_[blockI][faceI][nodeI];
				int global_node = local_2_global_nodes[blockI][local_node];
				face_2_nodes_connectivity.push_back(global_node);
			}


			std::vector<int> face_2_nodes_connectivity_non_sorted= face_2_nodes_connectivity;
			std::sort(face_2_nodes_connectivity.begin(),face_2_nodes_connectivity.end());
			std::vector<int>::iterator it;
			std::vector<int> face_intersection(face_2_nodes_connectivity.size()+common_nodes_vector.size());

			it = std::set_intersection( face_2_nodes_connectivity.begin(),face_2_nodes_connectivity.end(), common_nodes_vector.begin(), common_nodes_vector.end(), face_intersection.begin() );
			face_intersection.resize(it-face_intersection.begin());

			if(face_intersection.size()>=3)
			{
				common_faces_vector_temp[connexionI].push_back(face_2_nodes_connectivity_non_sorted);
				// commonCellsVector_[connexionI].push_back(local_2_global_elements[blockI][cell_local]);
				std::cout<<boundaryFaces_[blockI][faceI][0]<<"\t"<<boundaryFaces_[blockI][faceI][1]<<"\t"<<boundaryFaces_[blockI][faceI][2]<<"\t"<<boundaryFaces_[blockI][faceI][3]<<std::endl;
				
				std::cout<<face_2_nodes_connectivity[0]<<"\t"<<face_2_nodes_connectivity[1]<<"\t"<<face_2_nodes_connectivity[2]<<"\t"<<face_2_nodes_connectivity[3]<<std::endl;

				std::cout << "======================"<<common_faces_vector_temp[connexionI].size()<<"========================" << '\n';

			}
		}


		for(int faceI = 0;faceI<common_faces_vector_temp[connexionI].size();faceI++)
		{
			int node0 = common_faces_vector_temp[connexionI][faceI][0];
			int node1 = common_faces_vector_temp[connexionI][faceI][1];
			int node2 = common_faces_vector_temp[connexionI][faceI][2];

			std::vector<int> node_2_cells_connectivity_node0 = global_node_2_cells_connectivity[0][node0];
			std::vector<int> node_2_cells_connectivity_node1 = global_node_2_cells_connectivity[0][node1];
			std::vector<int> node_2_cells_connectivity_node2 = global_node_2_cells_connectivity[0][node2];

			std::sort(node_2_cells_connectivity_node0.begin(),node_2_cells_connectivity_node0.end());
			std::sort(node_2_cells_connectivity_node1.begin(),node_2_cells_connectivity_node1.end());
			std::sort(node_2_cells_connectivity_node2.begin(),node_2_cells_connectivity_node2.end());

			std::vector<int>::iterator first_it;
			std::vector<int> first_cells_intersection(node_2_cells_connectivity_node0.size()+node_2_cells_connectivity_node1.size());

			first_it = std::set_intersection( node_2_cells_connectivity_node0.begin(), node_2_cells_connectivity_node0.end(), node_2_cells_connectivity_node1.begin(), node_2_cells_connectivity_node1.end(), first_cells_intersection.begin() );
			first_cells_intersection.resize(first_it-first_cells_intersection.begin());

			if(first_cells_intersection.size()>0)
			{
				std::vector<int>::iterator second_it;
				std::vector<int> second_cells_intersection(node_2_cells_connectivity_node2.size()+first_cells_intersection.size());
				second_it = std::set_intersection( node_2_cells_connectivity_node2.begin(), node_2_cells_connectivity_node2.end(), first_cells_intersection.begin(), first_cells_intersection.end(), second_cells_intersection.begin() );

				second_cells_intersection.resize(second_it - second_cells_intersection.begin());

				if(second_cells_intersection.size()==2)
				{
					int cell0 = second_cells_intersection[0];
					int cell1 = second_cells_intersection[1];

					int block_cell0 = global_2_local_elements[cell0][1];
					int block_cell1 = global_2_local_elements[cell1][1];

					if((block_cell0==first_block_id || block_cell0 == second_block_id)&&(block_cell1==first_block_id || block_cell1 == second_block_id))
					{
						std::cout<<"ok \n";
						commonFacesVector_[connexionI].push_back(common_faces_vector_temp[connexionI][faceI]);
						commonCellsVector_[connexionI].push_back({cell0,cell1});			
					}


				}
			}

		}
	}
}

void ReconstructFaces::RemovePhysicalBoundaries(int** local_2_global_nodes, MetisBoundary* single_block_metisboundary)
{

	for(int blockI =0;blockI<n_blocks_;blockI++)
	{
		for(int faceI =0 ; faceI<boundaryFaces_[blockI].size(); faceI++)
		{
			std::vector<int> face_2_nodes_connectivity_local = boundaryFaces_[blockI][faceI];

			std::vector<int> face_2_nodes_connectivity_global;

			for(int nodeI=0;nodeI<boundaryFaces_[blockI][faceI].size();nodeI++)
			{
				int local_node = boundaryFaces_[blockI][faceI][nodeI];
				int global_node = local_2_global_nodes[blockI][local_node];
				face_2_nodes_connectivity_global.push_back(global_node);
			}
			std::vector<int> face_2_nodes_connectivity_global_non_sorted = face_2_nodes_connectivity_global;

			std::sort(face_2_nodes_connectivity_global.begin(),face_2_nodes_connectivity_global.end());
			
			bool found_face =false;

			for(int boundaryI = 0; boundaryI < single_block_metisboundary->nBoundaries_;boundaryI++)
			{
				for(int boundary_faceI =0; boundary_faceI < single_block_metisboundary->boundaryNelements_[boundaryI]; boundary_faceI++)
				{
					std::cout<<"n_elements_in_boundary = "<<single_block_metisboundary->boundaryNelements_[boundaryI]<<std::endl;

					std::cout<<"boundary_faceI = "<<boundary_faceI<<std::endl;

					std::vector<int> boundary_face_2_nodes_connectivity = single_block_metisboundary->boundaryConnectivity_[boundaryI][boundary_faceI];

					std::sort( boundary_face_2_nodes_connectivity.begin(), boundary_face_2_nodes_connectivity .end());
					std::vector<int>::iterator it;
					std::vector<int> faces_intersection(boundary_face_2_nodes_connectivity.size()+face_2_nodes_connectivity_global.size());

					it = std::set_intersection(boundary_face_2_nodes_connectivity.begin(), boundary_face_2_nodes_connectivity.end(),face_2_nodes_connectivity_global.begin(), face_2_nodes_connectivity_global.end(),faces_intersection.begin());
					faces_intersection.resize(it-faces_intersection.begin());

					if(faces_intersection.size()>=3)
					{
						boundaryFaces_[blockI][faceI].clear();
						found_face = true;
						break;
					}

				}
				if(found_face)
				{
					break;
				}
			}
		}
	}
}

