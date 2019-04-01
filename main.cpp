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
	cout << "========================STARTING PROGRAM========================" << endl;

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

if (argc != 4)
{
	std::cout << "Usage: ./metis <single block mesh file> <Number of partitions> <Output mesh file name>\n";
	return 0;
}

  //Metis' routine

  //Input arguments
std::string meshFile = argv[1];
int nPart = atoi(argv[2]);
std::string outputMeshFile = argv[3];


MetisMesh reader;
reader.ReadSingleBlockMesh(meshFile);

int n_blocks =nPart;
int* global_n_elements = reader.getNElements_();


MetisMesh* newMesh = reader.Partition(nPart);

ReconstructFaces reconstruct_faces = ReconstructFaces(n_blocks, global_n_elements[0]);
std::vector<int> common_nodes_vector;

int** node_array = newMesh ->getLocal2GlobalNodes_();
int* n_nodes_in_block = newMesh ->getNNodes_();

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
cout<<reconstruct_faces.connexionVector_.size()<<endl;
for(int i = 0; i<reconstruct_faces.connexionVector_.size();i++)
{
	for(int j = 0; j<reconstruct_faces.connexionVector_[i].size();j++)
	{
		cout<<reconstruct_faces.connexionVector_[i][j]<<"\t";
	}
	std::cout << '\n';
}

newMesh->WriteMesh(outputMeshFile);
newMesh->WriteOutputTecplot("outputMeshFile.dat");

}
