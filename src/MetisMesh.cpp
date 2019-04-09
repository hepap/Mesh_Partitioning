#include "MetisMesh.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <metis.h>
#include <map>
#include <string>
#include <algorithm>
#include "omp.h"

using namespace std;

/*==============HELENE=================*/
int* MetisMesh::getNNodes_()
{
  return nNodes_;
}
int* MetisMesh::getNElements_()
{
  return nElements_;
}
int** MetisMesh::getlocal2GlobalNodes_()
{
  return local2GlobalNodes_;
}
int** MetisMesh::getlocal2GlobalElements_()
{
  return local2GlobalElements_;
}
std::vector<int>* MetisMesh::getGlobal2LocalNodes_()
{
  return global2LocalNodes_;
}
std::vector<int>* MetisMesh::getGlobal2LocalElements_()
{
  return global2LocalElements_;
}
std::vector<int>** MetisMesh::getConnectivity_()
{
  return connectivity_;
}
std::vector<int>** MetisMesh::getNode2Cells_()
{
  return node2Cells_;
}

/*=====================================*/


/*==============Isabelle=================*/
MetisBoundary* MetisMesh::GetMetisBoundary_()
{
  return metisBoundary_;
}
int* MetisMesh::getElementBlock_()
{
  return elementBlock_;
}
/*=======================================*/


int findNodeIndex(std::vector<int> &list, int node2find)
{
    int count(list.size());

    for (int i = 0; i < count; i++)
    {
        if (node2find == list[i])
        {
            return i+1;
        }
    }

    list.push_back(node2find);
    return list.size();
}

MetisMesh::MetisMesh()
    : nElements_(nullptr), nNodes_(nullptr), elementNbrNodes_(nullptr), elementType_(nullptr),
     local2GlobalElements_(nullptr),local2GlobalNodes_(nullptr),global2LocalNodes_(nullptr),
     node2Cells_(nullptr), cell2GlobalNodes_(nullptr), nTotalNode_(0), nBlock_(0),
     x_(nullptr), y_(nullptr), z_(nullptr), connectivity_(nullptr), elementBlock_(nullptr),
     localBoundary_()



{
    std::cout << "constructing MetisMesh..." << endl;
}

MetisMesh::~MetisMesh()
{
    std::cout << "destructing MetisMesh..." << endl;

    if (nElements_ != nullptr)
        delete[] nElements_;
    if (nNodes_ != nullptr)
        delete[] nNodes_;

    nElements_ = nullptr;
    nNodes_ = nullptr;


    for (int blockI = 0; blockI < nBlock_; blockI++)
    {
        if (x_[blockI] != nullptr)
            delete[] x_[blockI];
        if (y_[blockI] != nullptr)
            delete[] y_[blockI];
        if (z_[blockI] != nullptr)
            delete[] z_[blockI];


        if (local2GlobalElements_[blockI] != nullptr)
            delete[] local2GlobalElements_[blockI];
        if (connectivity_[blockI] != nullptr)
            delete[] connectivity_[blockI];


        connectivity_[blockI] = nullptr;
        /*==============HELENE=================*/
        // if (node2Cells_[blockI] != nullptr)
        //     delete[] node2Cells_[blockI];
        //
        //
        // node2Cells_[blockI] = nullptr;
        /*=====================================*/
    }

    std::cout << "destructor ok 2" << endl;

    if (x_ != nullptr)
        delete[] x_;
    if (y_ != nullptr)
        delete[] y_;
    if (z_ != nullptr)
        delete[] z_;

    if (elementType_ != nullptr)
        delete[] elementType_;
     if (elementNbrNodes_ != nullptr)
        delete[] elementNbrNodes_;

    std::cout << "destructor ok 3" << endl;
    if (local2GlobalElements_ != nullptr)
        delete[] local2GlobalElements_;

    if (connectivity_ != nullptr)
        delete[] connectivity_;

    /*==============HELENE=================*/
    // if (node2Cells_ != nullptr)
    //     delete[] node2Cells_;
    /*=====================================*/


    x_ = nullptr;
    y_ = nullptr;
    z_ = nullptr;

    local2GlobalElements_ = nullptr;
    elementType_ = nullptr;
    elementNbrNodes_ = nullptr;
    connectivity_ = nullptr;

    /*==============HELENE=================*/
    // node2Cells_ = nullptr;
    /*=====================================*/


    std::cout << "destructor ok 4" << endl;
    std::cout << "destruct succeed" << endl;
}

void MetisMesh::Init(int nBlock, int *nElements, int *nNodes)

{
    std::cout << "starting Init..." << endl;
    nBlock_ = nBlock;
    nElements_ = new int[nBlock];
    nNodes_ = new int[nBlock];


    /*==============HELENE=================*/
    elementType_ = new int [nBlock];
    elementNbrNodes_ = new int [nBlock];
    node2Cells_ = new std::vector<int> *[nBlock];
    /*=====================================*/



    x_ = new double *[nBlock];
    y_ = new double *[nBlock];
    z_ = new double *[nBlock];

    connectivity_ = new std::vector<int> *[nBlock];



    for (int i = 0; i < nBlock; i++)
    {
        nElements_[i] = nElements[i];
        nNodes_[i] = nNodes[i];
        x_[i] = new double[nNodes[i]];
        y_[i] = new double[nNodes[i]];
        z_[i] = new double[nNodes[i]];
       /*  elementType_[i] = new int [nElements[i]];
        elementNbrNodes_[i] = new int [nElements[i]]; */

        connectivity_[i] = new std::vector<int>[nElements_[i]];
        /*==============HELENE=================*/
        node2Cells_[i] = new std::vector<int>[nNodes_[i]];
        /*=====================================*/

    }

    // Only want to initialize after partitioning
    if (nBlock != 1) {
        local2GlobalElements_ = new int *[nBlock];

        for (int i = 0; i < nBlock; i++) {
            local2GlobalElements_[i] = new int [nElements[i]];
        }
        std::cout << "Initialisation done global" << endl;
    }


}

// Initialize static int attribute nDimensions_
int MetisMesh::nDimensions_ = 0;

void MetisMesh::ReadSingleBlockMesh(std::string fileName)
{
    std::cout << "Start Reading Function " << endl;
    ifstream myfile(fileName);
    string line;

    if (myfile.is_open())
    {
        std::cout << "Reading...................... " << fileName << endl;

        int nNodes(0);
        int nElements(0);
        int nBoundaries(0);
        int nDimensions(0);
        int nBlock(1);
        char str_temp[100];

        // ------- Prepass 1 : -------- //
        getline(myfile, line);

        // Stocke le nombre de dimensions du maillage dans une variable statique
        sscanf(line.c_str(), "NDIME=%d", &nDimensions);
        std::cout << "Nombre de dimensions = " << nDimensions << endl;
        nDimensions_ = nDimensions;
        getline(myfile, line);

        // Stocke le nombre d'éléments du maillage
        sscanf(line.c_str(), "NELEM=%d", &nElements);
        std::cout << "Nombre d'elements= " << nElements << endl;

        for (int i = 0; i < nElements + 1; i++)
        {
            getline(myfile, line);
        }

        // Stocke le nombre de noeuds du maillage
        sscanf(line.c_str(), "NPOIN=%d", &nNodes);
        std::cout << "Nombre de noeuds= " << nNodes << endl;

        for (int i = 0; i < nNodes + 1; i++)
        {
            getline(myfile, line);
        }
        std::cout << line << endl;
        // Stocke le nombre de conditions frontieres
        sscanf(line.c_str(), "NMARK=%d", &nBoundaries);
        std::cout << "Nombre de conditions frontières= " << nBoundaries << endl;

        // Création d'un object MetisBoundary pour stocker les variables
        metisBoundary_ = new MetisBoundary(nBoundaries);


        for (int i = 0; i < nBoundaries; i++) {


            getline(myfile, line);
            // Finding the first equal sign in the string
            // This is because we want only its right side MARKER_TAG= AIRFOIL
            size_t pos_equal_sign = line.find('=');
            // The string tag is whatever is at the right of the equal sign.
            string str_tag = line.substr(pos_equal_sign+1);
            // The following is the equivalent to a left trim (strip)
            while (str_tag[0] == ' ')
            {
                str_tag = str_tag.substr(1);
            }

            metisBoundary_->boundaryNames_[i] = str_tag;

            std::cout << "le nom de la frontiere " << i+1 << " est : " << str_tag << endl;


            getline(myfile, line);

            int boundaryNelements = 0;

            sscanf(line.c_str(), "MARKER_ELEMS= %d", &boundaryNelements);

            // Objets a deleter
            int* elementTypeBoundary = new int[boundaryNelements];
            metisBoundary_->boundaryNelements_[i] = boundaryNelements;
            metisBoundary_->boundaryElementNbrNodes_[i] = new int[boundaryNelements];
            metisBoundary_->boundaryElementType_[i] = new int[boundaryNelements];

            for (int j = 0; j < boundaryNelements; j++) {

                getline(myfile, line);
                sscanf(line.c_str(), "%d %s", &elementTypeBoundary[j], str_temp);

                metisBoundary_->boundaryElementType_[i][j] = elementTypeBoundary[j];
                metisBoundary_->boundaryElementNbrNodes_[i][j] = NumberOfNodes(elementTypeBoundary[j]);
            }
        }

        // Fonctions d'initialisation :
        Init(nBlock, &nElements, &nNodes);
        metisBoundary_->InitBoundary(metisBoundary_->boundaryNelements_, nBoundaries);

        // Prepass 2 : retour au debut du fichier
        myfile.seekg(0, myfile.beg);

        getline(myfile, line);
        getline(myfile, line);
        getline(myfile, line);

        // Enregistrement du premier entier elementType
        nTotalNode_ = 0;
        elementType_ = new int[nElements];
        elementNbrNodes_ = new int[nElements];

        for (int i = 0; i < nElements; i++)
        {
            sscanf(line.c_str(), "%d %s", &elementType_[i], str_temp);
            elementNbrNodes_[i] = NumberOfNodes(elementType_[i]);
            nTotalNode_ += elementNbrNodes_[i];
            getline(myfile, line);

        }

        std::cout << "nTotalNode_ = " << nTotalNode_ << endl;
        std::cout << "nNodes = " << nNodes << endl;

        if (nDimensions == 3)
        {
            for (int nodeI = 0; nodeI < nNodes; nodeI++)
            {
                myfile >> x_[0][nodeI] >> y_[0][nodeI] >> z_[0][nodeI];
                getline(myfile, line);
            }

        }
        else
        {
            for (int nodeI = 0; nodeI < nNodes; nodeI++)
            {
                myfile >> x_[0][nodeI] >> y_[0][nodeI];
                getline(myfile, line);
                z_[0][nodeI] = 0.0;
            }
        }

        std::cout << "Noeuds enregistres " << endl;
        getline(myfile, line);
        getline(myfile, line);
        getline(myfile, line);

        int node;
        int temp;

        for (int i = 0; i < nBoundaries; i ++) {

            for (int j = 0; j < metisBoundary_->boundaryNelements_[i]; j++) {
                myfile >> temp;
                for (int k = 0; k < metisBoundary_->boundaryElementNbrNodes_[i][j]; k++)
                {

                    myfile >> node;
                    metisBoundary_->boundaryConnectivity_[i][j].push_back(node);
                }

                getline(myfile, line);
            }

            if (i != nBoundaries - 1)
            {
                getline(myfile, line);
                getline(myfile, line);
            }
        }


        std::cout << "Frontieres enregistrees " << endl;


        // getline and tokenize node in boundaries
        myfile.seekg(0, myfile.beg);
        //debut du fichier remplacer par une recher mot cle ^
        getline(myfile, line);
        getline(myfile, line);


        for (int i = 0; i < nElements; i++)
        {
            myfile >> temp;

            for (int numberOfNodes = 0; numberOfNodes < elementNbrNodes_[i]; numberOfNodes++)
            {

                int node;
                myfile >> node;
                connectivity_[0][i].push_back(node);
                /*==============HELENE=================*/
                node2Cells_[0][node].push_back(i);
                /*=====================================*/
            }

            getline(myfile, line);
        }

        myfile.close();
        std::cout << "closing... " << fileName << endl;

    }
    else
        std::cout << "could not open " << fileName << endl;

}



void MetisMesh::WriteMesh(std::string fileName, ReconstructFaces* reconstruct_faces)
{

    std::vector<string> filesName;

    for (int blockI = 0; blockI < nBlock_; blockI++)
    {
      int connexion_count =0;
      for(int connexionI=0;connexionI<reconstruct_faces->connexionVector_.size();connexionI++)
      {
        if ((reconstruct_faces->connexionVector_[connexionI][0]==blockI)||(reconstruct_faces->connexionVector_[connexionI][1]==blockI))
        {
          connexion_count +=1;
        }

      }



        // creer dynamiquement le nom des fichiers partition a partir de linput donne par le user
        std::string name = fileName + std::to_string(blockI) + ".su2"; // C++11 for std::to_string
        // store name in a vector
        filesName.push_back(name);

        FILE *fid = fopen(name.c_str(), "w");
        std::cout << "file open... " << name << endl;

        int nNodes = nNodes_[blockI];
        int nElements = nElements_[blockI];


        // fprintf(fid, "Block= %d\n", blockI);
        fprintf(fid, "NDIME= %d\n", nDimensions_);
        fprintf(fid, "\n");

        std::cout << "access connectivity elements block " << blockI << endl;
        fprintf(fid, "NPOIN= %d\n", nNodes);
        for (int nodeI = 0; nodeI < nNodes; nodeI++)
        {
            fprintf(fid, "%.12e %.12e %.12e\n", x_[blockI][nodeI], y_[blockI][nodeI], z_[blockI][nodeI]);
        }
        fprintf(fid, "\n");

        fprintf(fid, "NELEM= %d\n", nElements);

        for (int elementI = 0; elementI < nElements; elementI++)
        {

            int elementGlobal = local2GlobalElements_[blockI][elementI];

            //cout << "element local du block" << blockI << " element" << elementI <<  " = " << elementGlobal << endl;
            fprintf(fid, "%d ", elementType_[elementGlobal]);
            //cout << "elementType_[elementGlobal] " << elementType_[elementGlobal] << endl;
            //cout << "elementNbrNodes_[elementGlobal] " << elementNbrNodes_[elementGlobal] << endl;


           for (int j = 0; j < elementNbrNodes_[elementGlobal]; j++)
            {
                fprintf(fid, "%d ", connectivity_[blockI][elementI][j]);
            }

            fprintf(fid, "\n");
        }



        // Ajoute des frontiere physique
        std::cout << " ----- Ajout des frontieres pour le block " << blockI << endl;
        fprintf(fid, "\n");

        fprintf(fid, "NMARK= %d\n", metisBoundary_->nBoundaries_+connexion_count);
        std::cout << "Nombre de frontieres = " << metisBoundary_->nBoundaries_ << endl;

        for (int boundaryI = 0; boundaryI < metisBoundary_->nBoundaries_; boundaryI++)
        {

            int markerElems = (*localBoundary_)[make_pair(boundaryI, blockI)].size();
            std::cout << "markerElems " << markerElems << endl;

            fprintf(fid, "MARKER_TAG= %s\n", metisBoundary_->boundaryNames_[boundaryI].c_str());
            fprintf(fid, "MARKER_ELEMS= %d\n", markerElems);

            for (int elementI = 0; elementI < markerElems; elementI++) {


                int NELEM = (*localBoundary_)[make_pair(boundaryI, blockI)][elementI].size();
                int su2Id;

                switch(NELEM) {
                    case 4 : su2Id = 9;
                        break;
                    case 3 : su2Id = 5;
                        break;
                    case 2 : su2Id = 2;
                        break;
                }

                fprintf(fid, "%d ", su2Id);
                for (int nodeI = 0; nodeI < NELEM; nodeI ++) {
                    fprintf(fid, "%d ", (*localBoundary_)[make_pair(boundaryI, blockI)][elementI][nodeI]);
                   // std::cout << (*localBoundary_)[make_pair(boundaryI, blockI)][elementI][nodeI] << " ";
                }

                //std::cout << endl;
                fprintf(fid, "\n");
            }

        }

        for(int connexionI=0;connexionI<reconstruct_faces->connexionVector_.size();connexionI++)
        {
          if ((reconstruct_faces->connexionVector_[connexionI][0]==blockI)||(reconstruct_faces->connexionVector_[connexionI][1]==blockI))
          {
            fprintf(fid, "MARKER_TAG= CONNEXION\n");
            fprintf(fid, "MARKER_ELEMS= %d\n", reconstruct_faces->commonFacesVector_[connexionI].size());

            for(int faceI=0;faceI<reconstruct_faces->commonFacesVector_[connexionI].size(); faceI++)
            {
              std::vector<int> global_face_2_nodes_connectivity = reconstruct_faces->commonFacesVector_[connexionI][faceI];
              int face_type =0;
              if(global_face_2_nodes_connectivity.size()==3)
              {
                face_type = 5;
              }
              else if(global_face_2_nodes_connectivity.size()==4)
              {
                face_type = 9;
              }
              fprintf(fid, "%d", face_type);
              for(int nodeI =0 ; nodeI<global_face_2_nodes_connectivity.size();nodeI++)
              {
                int global_node = global_face_2_nodes_connectivity[nodeI];
                int local_node = ReturnLocalNode(global_node, blockI, global2LocalNodes_);
                fprintf(fid, " %d", local_node);
              }
              fprintf(fid, "\n");
            }
          }


        }

        fclose(fid);
        std::cout << name << " output file closed ..." << endl;
    }


    // store name in attribute so we can print them later in topology gile
    filesName_ = filesName;
}

int MetisMesh::NumberOfNodes(int elementType)
{

    int numberOfNodes(0);

    // Line
    if (elementType == 3)
    {

        numberOfNodes = 2;
    }

    // Triangle
    else if (elementType == 5)
    {

        numberOfNodes = 3;
    }

    // Quadrilateral
    else if (elementType == 9)
    {

        numberOfNodes = 4;
    }

    // Tetrahedral
    else if (elementType == 10)
    {

        numberOfNodes = 4;
    }

    // Hexahedral
    else if (elementType == 12)
    {

        numberOfNodes = 8;
    }

    // Wedge
    else if (elementType == 13)
    {

        numberOfNodes = 6;
    }

    // Pyramid
    else if (elementType == 14)
    {

        numberOfNodes = 5;
    }

    return numberOfNodes;
}

MetisMesh* MetisMesh::Partition(int nPart)
{
    std::cout << "fonction Partition start " << endl;

    if (nBlock_ > 1)
    {
        throw std::runtime_error("Cannot partition a multiblock mesh!");
    }

    // Vecteurs d'entrée pour Metis
    int eptr[nElements_[0] + 1];
    int *eind = new int[nTotalNode_];

    // Converting connectivity into METIS data structure See Metis reference doc
    eptr[0] = 0;
    int count = 0;
    int count1 = 0;

    for (int i = 0; i < nElements_[0]; i++)
    {

        count += elementNbrNodes_[i];
        eptr[i+1] = count; // i + 1, car eptr[0] doit commencer a 0
    }


    for (int i = 0; i < nElements_[0]; i++)
    {
        for (int j = 0; j < elementNbrNodes_[i]; j++)
        {
            eind[count1] = connectivity_[0][i][j];
            count1++;

        }
    }

    std::cout << "fin connectivite Metis" << endl;

    int ncommon(3);
    int objval;
    int epart[nElements_[0]];
    int npart[nNodes_[0]];

    int success = METIS_PartMeshDual(&nElements_[0], &nNodes_[0], &eptr[0], eind, NULL, NULL,
                                     &ncommon, &nPart, NULL, NULL, &objval,
                                     &epart[0], &npart[0]);

    std::cout << "Partition Success: " << success << std::endl;


    // for(int i=0;i<16;i++)
    // {
    //   epart[i] = 0;
    //   epart[i+16] = 1;
    //   epart[i+32] = 2;
    //   epart[i+48] = 3;

    // }


    std::vector<int> elementsPerBlock[nPart];
    std::vector<int> elementNbrNodesPerBlock[nPart];


    int* elementBlock_ = new int[nElements_[0]];

    for (int i = 0; i < nElements_[0]; i++)
    {
        elementBlock_[i] = epart[i];
        //std::cout << "element" << i << " = " << elementBlock_[i] << endl;
        int blockId = epart[i];

        elementsPerBlock[blockId].push_back(i);
        elementNbrNodesPerBlock[blockId].push_back(elementNbrNodes_[i]);
    }


    /*==============HELENE=================
    std::vector<int> nodesPerBlock[nPart];
    for (int i = 0; i < nNodes_[0]; i++)
    {
        int blockId = npart[i];
        nodesPerBlock[blockId].push_back(i);
        elementNbrNodesPerBlock[blockId].push_back(elementNbrNodes_[i]);
    }
    =====================================*/

    int newNelements[nPart];

    for (int blockI = 0; blockI < nPart; blockI++) {
        newNelements[blockI] = elementsPerBlock[blockI].size();
    }

    std::cout << "Partition and new vector" << endl;

    std::vector<int> addedNode[nPart];
    std::vector<int> **newConnectivity;
    newConnectivity = new std::vector<int> *[nPart];

    /*==============HELENE=================*/
    // std::vector<int> **newConnectivity_node;
    // newConnectivity_node = new std::vector<int> *[nPart];
    /*=====================================*/

    std::cout << "nPart = " << nPart << endl;

    // #pragma omp parallel for num_threads(4)
    for (int blockI = 0; blockI < nPart; blockI++)
    {
        newConnectivity[blockI] = new std::vector<int>[newNelements[blockI]];


        std::cout << "building newConnectivity " << blockI << endl;
        std::cout << newNelements[blockI] << endl;

        for (int i = 0; i < newNelements[blockI]; i++)
        {
            int elementIblockI = elementsPerBlock[blockI][i];
            int vector_size = NumberOfNodes(elementType_[elementIblockI]);
            newConnectivity[blockI][i].resize(vector_size);
        }

        for (int i = 0; i < newNelements[blockI]; i++)
        {
            for (int j = 0; j < elementNbrNodesPerBlock[blockI][i]; j++)
            {
                int elementIblockI = elementsPerBlock[blockI][i];
                int n1 = connectivity_[0][elementIblockI][j];
                int newN1 = findNodeIndex(addedNode[blockI], n1);
                newConnectivity[blockI][i][j] = newN1-1;

            }
        }

    }
    std::cout << "new connectivity done" << endl;
    int newNnodes[nPart];
    /*==============HELENE=================*/
    int nodeCount4Global2LocalNodes = 0;
    /*=====================================*/
    // #pragma omp parallel for num_threads(4)
    for (int blockI = 0; blockI < nPart; blockI++)
    {
        newNnodes[blockI] = addedNode[blockI].size();

        /*==============HELENE=================*/
        nodeCount4Global2LocalNodes += newNnodes[blockI];
        /*=====================================*/
    }

    // Creation dune nouvelle instance de Metis Mesh
    MetisMesh *newMesh = new MetisMesh();
    newMesh->Init(nPart, newNelements, newNnodes);

    newMesh->SetConnectivity(newConnectivity);

    std::cout << "newMesh ok" << endl;
    newMesh->node2Cells_ = new std::vector<int> *[nPart];

    for (int blockI = 0; blockI < nPart; blockI++)
    {
      newMesh->node2Cells_[blockI] = new std::vector<int>[newMesh->nNodes_[blockI]];
      for (int elementI = 0; elementI < newMesh->nElements_[blockI]; elementI++)
      {
        for(int nodeI=0; nodeI < newMesh->connectivity_[blockI][elementI].size();nodeI++)
        {
          int node = newMesh->connectivity_[blockI][elementI][nodeI];

          newMesh->node2Cells_[blockI][node].push_back(elementI);
        }
      }

    }






    newMesh->global2LocalElements_ = new std::vector<int>[nElements_[0]];

    // #pragma omp parallel for num_threads(4)
    for (int blockI = 0; blockI < nPart; blockI++)
    {
        newNelements[blockI] = elementsPerBlock[blockI].size();

        for (int i = 0; i < newNelements[blockI]; i++)
        {
            newMesh->local2GlobalElements_[blockI][i] = elementsPerBlock[blockI][i];

            newMesh->global2LocalElements_[elementsPerBlock[blockI][i]].push_back(i);
            newMesh->global2LocalElements_[elementsPerBlock[blockI][i]].push_back(blockI);

        }
    }

    std::cout << "openmp" << endl;
    /*==============HELENE=================*/
    newMesh->local2GlobalNodes_ = new int*[nPart];
    newMesh->global2LocalNodes_ = new std::vector<int>[nNodes_[0]];

    // for (int i = 0; i < nNodes_[0]; i++)
    // {
    //     newMesh->global2LocalNodes_[i] = new std::vector<int>;
    // }
    /*=====================================*/

    std::cout << "local and global" << endl;
    for (int blockI = 0; blockI < nPart; blockI++)
    {

        /*==============HELENE=================*/
        int size = addedNode[blockI].size();
        newMesh->local2GlobalNodes_[blockI] = new int[size];

        /*=====================================*/

        for (int i = 0; i < size; i++)
        {
            int nodeId = addedNode[blockI][i];
            newMesh->x_[blockI][i] = x_[0][nodeId];
            newMesh->y_[blockI][i] = y_[0][nodeId];
            newMesh->z_[blockI][i] = z_[0][nodeId];

            /*==============HELENE=================*/
            newMesh->local2GlobalNodes_[blockI][i] = nodeId;
            newMesh->global2LocalNodes_[nodeId].push_back(i);
            newMesh->global2LocalNodes_[nodeId].push_back(blockI);

            /*=====================================*/
        }
    }





    std::cout << "addedNode ok" << endl;

    for (int blockI = 0; blockI < nPart; blockI++)
    {
        if (newConnectivity[blockI] != nullptr)
            delete[] newConnectivity[blockI];
        newConnectivity[blockI] = nullptr;
    }

    if (newConnectivity != nullptr)
        delete[] newConnectivity;
    newConnectivity = nullptr;

    if (eind != nullptr)
        delete[] eind;
    eind = nullptr;



    // TODO: Check if newMesh should (like in this case), have a pointer to local2global or create its own structure.
    // Pas le choix de copier les attributs dans le newMesh car cest sur lui que l<on call la fonction WriteMesh
    newMesh->elementType_ = elementType_;
    local2GlobalElements_ = newMesh->local2GlobalElements_;
    newMesh->elementNbrNodes_ = elementNbrNodes_;
    newMesh->elementBlock_ = elementBlock_;
    newMesh->metisBoundary_ = metisBoundary_;

    std::cout << "newMesh returned" << endl;

    return newMesh;

}

void MetisMesh::SetConnectivity(std::vector<int> **connectivity)
{
    std::cout << "set connectivity function called" << endl;
    for (int blockI = 0; blockI < nBlock_; blockI++)
        for (int elementI = 0; elementI < nElements_[blockI]; elementI++)
            connectivity_[blockI][elementI] = connectivity[blockI][elementI];
}

/* void print_vector(string name, vector<int> vec) {
    cout << "Vector " << name << ": " << endl;
    for (auto c : vec){
        cout << c << " ";
    }
    cout << endl;
} */

void MetisMesh::ComputePhysicalBoundaries(MetisBoundary* metisBoundary, std::vector<int>** globalNode2GlobalCells, std::vector<int>* global2LocalNodes)
{
   // vector <int>** globalNode2LocalNodes = this->getGlobal2LocalNodes_();

    //std::vector<int> **newBoundaryElementType;
    //newBoundaryElementType = new std::vector<int> *[nBlock_];
    std::vector<int> node1;
    std::vector<int> node2;
    int commonGlobalElement;
    std::map<pair<int, int>, vector<vector<int>>>* localBoundary = new map<pair<int, int>, vector<vector<int>>>;

    for (int boundaryI = 0; boundaryI < metisBoundary->nBoundaries_; boundaryI++)
    {

        std::vector<int> face2Block;
        int boundaryElementNbr = metisBoundary->boundaryNelements_[boundaryI];

        //newBoundaryElementType[boundaryI] = new std::vector<int> [boundaryElementNbr];
        // Map <BoundaryIndex, BlockNumber> to a vector<vector<int>>

        // accede a chacun des faces de la frontiere
        for (int i = 0; i < boundaryElementNbr; i++)
        {

            int sizeBoundary = metisBoundary->boundaryElementNbrNodes_[boundaryI][i];
            int firstNode = metisBoundary->boundaryConnectivity_[boundaryI][i][0];
            node1 = globalNode2GlobalCells[0][firstNode];
            std::sort(node1.begin(), node1.end());
            int node1size = node1.size();

            // accede aux noeuds formant la face a la frontiere
            for (int j = 1; j < sizeBoundary; j++)
            {

                // un seul noeud de lelement face
                int globalNode = metisBoundary->boundaryConnectivity_[boundaryI][i][j];

                node2 = globalNode2GlobalCells[0][globalNode];
                std::sort(node2.begin(), node2.end());
                std::vector<int>::iterator it;
                std::vector<int> vCommon(node2.size() + node1size);

                it = std::set_intersection(node1.begin(), node1.end(), node2.begin(), node2.end(), vCommon.begin() );
				        vCommon.resize(it - vCommon.begin());

                if (vCommon.size() == 1)
                {
                    commonGlobalElement = vCommon[0];
                    //cout << "Lelement commun pour la face " << i << " est : " << commonGlobalElement << endl;
                    // structure contenant pour chaque face, dans quel block elles sont
                    face2Block.push_back(elementBlock_[commonGlobalElement]);
                    // std::cout << "La face a la frontiere " << i << " appartient au block " << face2Block[i] << endl;

                    //cout << "face2Block[i] " << face2Block[i] << endl;
                    //cout << "Lelement commun est dans le block " << elementBlock_[commonGlobalElement] << endl;
                    break;
                }
                else
                {
                    node1 = vCommon;
                }
                //std::cout << "La face a la frontiere " << i << "appartient au block " << face2Block[i] << endl;
            }


            int localBlock = face2Block[i];
            vector<int> localBoundaryElement;

            for (int k = 0; k < sizeBoundary; k++)
            {

                int globalNode = metisBoundary->boundaryConnectivity_[boundaryI][i][k];

                int localNode = ReturnLocalNode(globalNode, localBlock, global2LocalNodes);
                localBoundaryElement.push_back(localNode);

            }

            (*localBoundary)[make_pair(boundaryI, localBlock)].push_back(localBoundaryElement);
        }

    }

    // Allocate as attribute
    localBoundary_ = localBoundary;
}

int MetisMesh::ReturnLocalNode(int globalNode, int localBlock, std::vector<int>* global2LocalNodes)
{

    int localNode = 0;
    int i = 1;
    int index = global2LocalNodes[globalNode][i];

    while (index != localBlock) {
        i = i + 2;
        index = global2LocalNodes[globalNode][i];
    }

    localNode = global2LocalNodes[globalNode][i-1];

    return localNode;

}


void MetisMesh::WriteTopology(std::string fileName, ReconstructFaces* reconstruct_faces)
{

    FILE *fid = fopen(fileName.c_str(), "w");
    std::cout << filesName_.size() << endl;

    std::cout << "Topology file opened.................. " << fileName << endl;
    //std::cout << newMesh->filesName_[0] << endl;
    fprintf(fid, "NBLOCK= %d\n", nBlock_);

    for (int blockI = 0; blockI < nBlock_; blockI++) {
        fprintf(fid, "%s\n", filesName_[blockI].c_str());
    }

    for (int blockI = 0; blockI < nBlock_; blockI++)
    {


        fprintf(fid, "Block= %d\n", blockI );

        int n_ghost_cells=0;
        for (int boundaryI = 0; boundaryI < metisBoundary_->nBoundaries_; boundaryI++)
        {

            n_ghost_cells += (*localBoundary_)[make_pair(boundaryI, blockI)].size();
        }

        n_ghost_cells+=nElements_[blockI];
        fprintf(fid, "NGhost= %d\n", n_ghost_cells );

        // Check how many connexions in block
        vector<int> connexions_in_block_idx;

        for (int j=0;j<reconstruct_faces->connexionVector_.size();j++)
        {
            if ((reconstruct_faces->connexionVector_[j][0]==blockI)||(reconstruct_faces->connexionVector_[j][1]==blockI))
            {
                connexions_in_block_idx.push_back(j);
            }
        }

        // Print n_connexions_in_block
        fprintf(fid, "Nconnexion= %d\n", connexions_in_block_idx.size() );


        for (int j=0;j<connexions_in_block_idx.size();j++)
        {
            int connexion_idx=connexions_in_block_idx[j];
            // Print Nmark for connexions
            if (reconstruct_faces->connexionVector_[connexion_idx][0]==blockI)
            {
                fprintf(fid, "Nmark= %d\n", reconstruct_faces->connexionVector_[connexion_idx][1] );

            }
            else
            {
                fprintf(fid, "Nmark= %d\n", reconstruct_faces->connexionVector_[connexion_idx][0] );
            }

            // Print Nelems for connexions
            fprintf(fid, "Nelems= %d\n", reconstruct_faces->commonCellsVector_[connexion_idx].size() );

            // Print each elem
            for (int k=0;k<reconstruct_faces->commonCellsVector_[connexion_idx].size();k++)
            {
                // FOR NOW PRINT 0 AND 1 LOCAL ALONG WITH IDX BLOCK
                int global2local_elem_0=global2LocalElements_[reconstruct_faces->commonCellsVector_[connexion_idx][k][0]][0];
                int global2local_block_0=global2LocalElements_[reconstruct_faces->commonCellsVector_[connexion_idx][k][0]][1];
                int global2local_elem_1=global2LocalElements_[reconstruct_faces->commonCellsVector_[connexion_idx][k][1]][0];
                int global2local_block_1=global2LocalElements_[reconstruct_faces->commonCellsVector_[connexion_idx][k][1]][1];

                if (global2local_block_0==blockI)
                {
                    fprintf(fid, "%d\n", global2local_elem_1);
                }
                else
                {
                    fprintf(fid, "%d\n", global2local_elem_0);
                }
            }
        }
    }

    fclose(fid);
    std::cout << fileName << "output file closed ..." << endl;
}



void MetisMesh::WriteOutputTecplot(std::string fileName, int** node_flag, int** cell_flag)
{
    // FILE *fid = fopen(fileName.c_str(), "w");
    std::ofstream fid;
    fid.open(fileName);
    std::cout << "file open... " << fileName << endl;

    fid << "TTILE = \"Vizualisation of the partitioned mesh\""<<endl;
    fid << "VARIABLES=\"X\",\"Y\",\"Z\",\"FLAGNODE\",\"FLAGCELL\"" << endl;


    for (int blockI = 0; blockI < nBlock_; blockI++)
    {
        int nNodes = nNodes_[blockI];
        int nElements = nElements_[blockI];
        fid << "ZONE T=\"element"<<blockI <<"\""<< endl;
        fid << "Nodes=" << nNodes << ", " << "Elements=" << nElements << ", " << "ZONETYPE=FETETRAHEDRON" << endl;
        fid << "DATAPACKING=BLOCK" << endl;
        fid << "VARLOCATION=([5]=CELLCENTERED)" << endl;

        for (int nodeI = 0; nodeI < nNodes; nodeI++)
        {
            fid <<  x_[blockI][nodeI] <<endl;
        }
        for (int nodeI = 0; nodeI < nNodes; nodeI++)
        {
            fid <<  y_[blockI][nodeI] <<endl;
        }
        for (int nodeI = 0; nodeI < nNodes; nodeI++)
        {
            fid <<  z_[blockI][nodeI] <<endl;
        }
        for (int nodeI = 0; nodeI < nNodes; nodeI++)
        {
            fid <<  node_flag[blockI][nodeI] <<endl;
        }
        for (int elementI = 0; elementI < nElements; elementI++)
        {
            fid <<  cell_flag[blockI][elementI] <<endl;
        }

        for (int elementI = 0; elementI < nElements; elementI++)
        {
            int elementGlobal = local2GlobalElements_[blockI][elementI];

          for (int j = 0; j < elementNbrNodes_[elementGlobal]; j++)
          {
            fid << connectivity_[blockI][elementI][j]+1<<"\t";
          }
            fid<< "\n";
        }
    }

    fid.close();
    std::cout << fileName << "outputTecplot file closed ..." << endl;
}
