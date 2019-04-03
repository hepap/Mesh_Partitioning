#ifndef MESH_PARTITIONING_HEAD_METISMESH_H
#define MESH_PARTITIONING_HEAD_METISMESH_H
#include <vector>
#include <string>
#include "MetisBoundary.h"
#include "ReconstructFaces.h"
#include "omp.h"

class MetisMesh
{
private:
    int* nElements_;
    int* nNodes_;

    /*==============HELENE=================*/
    int* elementNbrNodes_;
    int* elementType_;
    /*=====================================*/

    int** local2GlobalElements_;

    /*==============HELENE=================*/
    int** local2GlobalNodes_;
    std::vector<int>* global2LocalNodes_;
    std::vector<int>** node2Cells_;
    std::vector<int>** cell2GlobalNodes_;
    /*=====================================*/

    int nTotalNode_;
    int nBlock_;
     //A INIT ET DESTRUCT

    double** x_;
    double** y_;
    double** z_;

    MetisBoundary* metisBoundary_;
    std::vector<int>** connectivity_;


    std::vector<int>* global2LocalElements_;
    //std::vector<int>** connectivity_boundary;

public:
    MetisMesh();
    ~MetisMesh();

    friend class MetisBoundary;
    static int nDimensions_;
    vector<vector<int>>* ReturnFaces(int blockI, int element);


public:
    MetisBoundary* GetMetisBoundary_();
    int* GetEpart_();
    void Init(int nBlock, int* nElements, int* nNodes);
    void ReadSingleBlockMesh(std::string fileName);
    void WriteMesh(std::string fileName);
    MetisMesh* Partition(int nPart);
    void SetConnectivity(std::vector<int>** connectivity_);
    int NumberOfNodes(int elementType);
    int* getNNodes_();
    int* getNElements_();
    int** getLocal2GlobalNodes_();
    std::vector<int>* getGlobal2LocalNodes_();
    std::vector<int>* getGlobal2LocalElements_();
    std::vector<int>** getConnectivity_();
    std::vector<int>**  getNode2Cells_();
    std::vector<std::string> filesName_;
    // Fonction WriteTopology
    void WriteTopology(std::string filename, ReconstructFaces* reconstruct_faces);
    void ComputePhysicalBoundaries(MetisBoundary* metisBoundary, std::vector<int>** globalNode2globalCells);
    void WriteOutputTecplot(std::string fileName, int** node_flag,int** cell_flag);

};

#endif
