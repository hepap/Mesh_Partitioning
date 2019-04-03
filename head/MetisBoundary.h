#ifndef MESH_PARTITIONING_HEAD_METISBOUNDARY_H
#define MESH_PARTITIONING_HEAD_METISBOUNDARY_H
#include <string>
#include <vector>

using namespace std;

enum BoundaryType {
    WALL = 0,
    FAR_FIELD = 1,
    CONNECT = 2
};

class MetisBoundary
{

  


  static BoundaryType FindBoundaryTypeFromTagStr(string tag_str);

  public:
    MetisBoundary(int nBoundaries);
    //MetisBoundary(int nElements, int* elementType, int* elementNbrNodes, int **boundaryElements);
    ~MetisBoundary();
    void InitBoundary(int* nElements, int nBoundaries);
    
    // ====== Isabelle ----
    // Compute les frontiere physique et prend en entree un pointeur vers un objet MetisBoundary
    //void ComputePhysicalBoundaries(MetisBoundary* metisBoundary);
    string SetBoundaryName(string boundaryName);

    //vector<vector<int>> localNodes;


    //int** boundariesConnectivity_;


    int nBoundaries_;
    string* boundaryNames_;
    int** boundaryElements_;
    std::vector<int>** boundaryConnectivity_;
    //int *nNodes_;

    int* boundaryNelements_;
    int** boundaryElementType_;
    int** boundaryElementNbrNodes_;

    //int blockID_;
    //int** boundaryType_;


};

#endif
