#include "MetisMesh.h"
#include "ReconstructFaces.h"
#include <iostream>
#include <string>
#include <mpi.h>
#include <metis.h>
#include <vector>

using namespace std;

int main(int argc, char* argv[])
{
/* 	cout << "========================STARTING PROGRAM========================" << endl;
	cout << R"(
*        (       )   (            (          (          (       )     ) (       )
  (  `       )\ ) ( /(   )\ )   (     )\ )  *   ))\ )  *   ))\ ) ( /(  ( /( )\ ) ( /( (
  )\))(  (  (()/( )\()) (()/(   )\   (()/(` )  /(()/(` )  /(()/( )\()) )\()|()/( )\()))\ )
 ((_)()\ )\  /(_)|(_)\   /(_)|(((_)(  /(_))( )(_))(_))( )(_))(_)|(_)\ ((_)\ /(_)|(_)\(()/(
 (_()((_|(_)(_))  _((_) (_))  )\ _ )\(_)) (_(_()|_)) (_(_()|_))   ((_) _((_|_))  _((_)/(_))_
 |  \/  | __/ __|| || | | _ \ (_)_\(_) _ \|_   _|_ _||_   _|_ _| / _ \| \| |_ _|| \| (_)) __|
 | |\/| | _|\__ \| __ | |  _/  / _ \ |   /  | |  | |   | |  | | | (_) | .` || | | .` | | (_ |
 |_|  |_|___|___/|_||_| |_|   /_/ \_\|_|_\  |_| |___|  |_| |___| \___/|_|\_|___||_|\_|  \___|

	)" << endl;


 */

/*=============================TEST HELENE=============================
		int n_blocks = 2;
		int global_n_elements = 10;

ReconstructFaces reconstruct_faces = ReconstructFaces(n_blocks, global_n_elements);


std::vector<int> common_nodes_vector;
int first_block_id = 0;
int* first_block_node_array = new int[8] {0,1,2,3,4,5,8,9};
int n_nodes_in_first_block = 8;

int second_block_id =1;
int* second_block_node_array = new int[6] {1,2,3,4,6,7};
int n_nodes_in_second_block = 6;


common_nodes_vector = reconstruct_faces.CompareArraysOfGlobalNodes(  first_block_id,  first_block_node_array,  n_nodes_in_first_block,  second_block_id,  second_block_node_array, n_nodes_in_second_block);

std::vector<int>** global_cells_vector =

FindElementsInConnexion(std::vector<int> common_nodes_vector, std::vector<int>** global_cells_vector, std::vector<int>** global_nodes_vector, int block_id )


reconstruct_faces.~ReconstructFaces();
=============================TEST HELENE=============================*/

if (argc != 5)
{
	std::cout << "Usage: ./metis <single block mesh file> <Number of partitions> <Output mesh file name> <Output topology file name>\n";
	return 0;
}

  //Metis' routine

  //Input arguments
std::string meshFile = argv[1];
int nPart = atoi(argv[2]);
std::string outputMeshFile = argv[3];
std::string outputTopology = argv[4];

MetisMesh reader;
reader.ReadSingleBlockMesh(meshFile);
cout << "reading done" << endl;
int n_blocks =nPart;
int* global_n_elements = reader.getNElements_();


MetisMesh* newMesh = reader.Partition(nPart);

ReconstructFaces reconstruct_faces = ReconstructFaces(n_blocks, global_n_elements[0]);
std::vector<int> common_nodes_vector;

int** node_array = newMesh ->getLocal2GlobalNodes_();
int* n_nodes_in_block = newMesh ->getNNodes_();
int* n_elements_in_block = newMesh ->getNElements_();


for(int k=0; k<n_blocks; k++)
{
	for(int l=0; l<n_blocks;l++)
	{
		int first_block_id = k;
		int second_block_id = l;
		if(reconstruct_faces.block_array_4_comparaison_[k][l]!=1)
		{
			reconstruct_faces.CompareArraysOfGlobalNodes(  first_block_id,  node_array[k],  n_nodes_in_block[k],  second_block_id,  node_array[l], n_nodes_in_block[l]);
		}
	}
}



// std::vector<std::vector<int>> connexionVector = reconstruct_faces.connexionVector_;
// cout<<connexionVector.size()<<endl;
// for(int i= 0;i<connexionVector.size();i++)
// {
// 	for(int j=0;j<connexionVector[i].size();j++)
// 	{
// 		cout<<connexionVector[i][j]<<"\t";
// 	}
// 	cout<<"\n";
// }

std::vector<int>** globalCell2GlobalNodes = reader.getConnectivity_();
std::vector<int>** globalNode2GlobalCells = reader.getNode2Cells_();
int* elementBlock = reader.getElementBlock_();
reconstruct_faces.FindElementsInConnexion(globalCell2GlobalNodes,globalNode2GlobalCells);

cout<<"out of my face hehe"<<endl;
std::vector<std::vector<int>>* face2nodes = reconstruct_faces.commonFacesVector_;
cout<<"Is this happening?!"<<endl;
// cout<<face2nodes[0].size()<<endl;
// for(int i= 0;i<face2nodes[0].size();i++)
// {
// 	for(int j=0;j<face2nodes[0][i].size();j++)
// 	{
// 		cout<<face2nodes[0][i][j]<<"\t";
// 	}
// 	cout<<"\n";
// }
/*=============================TEST HELENE=============================*/
int** node_flag = new int*[n_blocks];
int** cell_flag = new int*[n_blocks];

std::vector<int>* global2LocalNodes =newMesh->getGlobal2LocalNodes_();
std::vector<int>* global2LocalElements =newMesh->getGlobal2LocalElements_();


for(int i=0;i<n_blocks;i++)
{
	node_flag[i] = new int[n_nodes_in_block[i]]();
	cell_flag[i] = new int[n_elements_in_block[i]]();
}

int n_connexions = reconstruct_faces.connexionVector_.size();
//
for(int i = 0; i<n_connexions;i++)
{
	for(int j = 0; j<reconstruct_faces.commonFacesVector_[i].size();j++)
	{
		int globalcell0 = reconstruct_faces.commonCellsVector_[i][j][0];
		int cell0 = global2LocalElements[globalcell0][0];
		int cellblock0 = global2LocalElements[globalcell0][1];
		cell_flag[cellblock0][cell0] = 2;
		int globalcell1 = reconstruct_faces.commonCellsVector_[i][j][1];
		int cell1 = global2LocalElements[globalcell1][0];
		int cellblock1 = global2LocalElements[globalcell1][1];
		cell_flag[cellblock1][cell1] = 2;
	/* 	cout<<reconstruct_faces.commonCellsVector_[i][j][0]<<"\t";
		cout<<cellblock0<<"\t";
		cout<<reconstruct_faces.commonCellsVector_[i][j][1]<<"\t";
		cout<<cellblock1<<"\t"; */

		for(int k = 0; k<reconstruct_faces.commonFacesVector_[i][j].size();k++)
		{
			//cout<<reconstruct_faces.commonFacesVector_[i][j][k]<<"\t";
			int global_node = reconstruct_faces.commonFacesVector_[i][j][k];
			if(global2LocalNodes[global_node].size()==4)
			{
				int local_node0 = global2LocalNodes[global_node][0];
				int block0 = global2LocalNodes[global_node][1];
				int local_node1 = global2LocalNodes[global_node][2];
				int block1 = global2LocalNodes[global_node][3];
				node_flag[block0][local_node0] = 1;
				node_flag[block1][local_node1] = 1;
			}
			else if(global2LocalNodes[global_node].size()==6)
			{
				int local_node0 = global2LocalNodes[global_node][0];
				int block0 = global2LocalNodes[global_node][1];
				int local_node1 = global2LocalNodes[global_node][2];
				int block1 = global2LocalNodes[global_node][3];
				int local_node2 = global2LocalNodes[global_node][4];
				int block2 = global2LocalNodes[global_node][5];
				node_flag[block0][local_node0] = 1;
				node_flag[block1][local_node1] = 1;
				node_flag[block2][local_node2] = 1;
			}
			else if(global2LocalNodes[global_node].size()==8)
			{
				int local_node0 = global2LocalNodes[global_node][0];
				int block0 = global2LocalNodes[global_node][1];
				int local_node1 = global2LocalNodes[global_node][2];
				int block1 = global2LocalNodes[global_node][3];
				int local_node2 = global2LocalNodes[global_node][4];
				int block2 = global2LocalNodes[global_node][5];
				int local_node3 = global2LocalNodes[global_node][6];
				int block3 = global2LocalNodes[global_node][7];
				node_flag[block0][local_node0] = 1;
				node_flag[block1][local_node1] = 1;
				node_flag[block2][local_node2] = 1;
				node_flag[block3][local_node3] = 1;
			}
			else if(global2LocalNodes[global_node].size()==10)
			{
				int local_node0 = global2LocalNodes[global_node][0];
				int block0 = global2LocalNodes[global_node][1];
				int local_node1 = global2LocalNodes[global_node][2];
				int block1 = global2LocalNodes[global_node][3];
				int local_node2 = global2LocalNodes[global_node][4];
				int block2 = global2LocalNodes[global_node][5];
				int local_node3 = global2LocalNodes[global_node][6];
				int block3 = global2LocalNodes[global_node][7];
				int local_node4 = global2LocalNodes[global_node][8];
				int block4 = global2LocalNodes[global_node][9];
				node_flag[block0][local_node0] = 1;
				node_flag[block1][local_node1] = 1;
				node_flag[block2][local_node2] = 1;
				node_flag[block3][local_node3] = 1;
				node_flag[block4][local_node4] = 1;

			}
		}
		//std::cout << endl;

	}

}
/*=============================TEST HELENE=============================*/

// ===== Isabelle =====

MetisBoundary* metisBoundary = reader.GetMetisBoundary_();

cout << "ComputeBoundaries " << endl;

newMesh->ComputePhysicalBoundaries(metisBoundary, globalNode2GlobalCells, global2LocalNodes);
// ====================

newMesh->WriteMesh(outputMeshFile);

// ===== Isabelle =====
newMesh->WriteOutputTecplot("outputMeshFile.dat",node_flag,cell_flag);
newMesh->WriteTopology(outputTopology, &reconstruct_faces);

}
