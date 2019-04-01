#include "MetisBoundary.h"
#include <iostream>
#include <string>

MetisBoundary::MetisBoundary(int nBoundaries) 
: nBoundaries_(nBoundaries), boundaryNames_(new string[nBoundaries]), boundaryConnectivity_(nullptr), boundaryNelements_(new int[nBoundaries]), 
boundaryElementType_(new int* [nBoundaries]), boundaryElementNbrNodes_(new int* [nBoundaries]) 
{
    cout << "constructing MetisBoundary........." << endl;
}


/* MetisBoundary::MetisBoundary(int nElements, int* elementType, int* elementNbrNodes, int **boundaryElements)
    : nElements_(nElements), elementType_(elementType),  
    elementNbrNodes_(elementNbrNodes),  boundaryElements_(boundaryElements) {}
          */

MetisBoundary::~MetisBoundary() {
    //if (nNodes_ != nullptr)
    //    delete[] nNodes_;
    //nNodes_ = nullptr;
}

 
void MetisBoundary::InitBoundary(int* nElements, int nBoundaries) {

  
    boundaryConnectivity_ = new std::vector<int> *[nBoundaries];

    for (int i = 0; i < nBoundaries; i++) {
        boundaryConnectivity_[i] = new std::vector<int>[nElements[i]];
    }

    

} 

// Arrive pas a etre callee depuis le main
/* void MetisBoundary::ComputePhysicalBoundaries(MetisBoundary* metisBoundary) 
{
    cout << "ALLLLOOO " << endl;

    for (int i = 0; i < metisBoundary->nBoundaries_; i++) {
        
        for (int j = 0; j < metisBoundary->boundaryNelements_[i]; j++) {
                
            for (int k = 0; k < metisBoundary->boundaryElementNbrNodes_[i][j]; k++)
            {
   
                cout << boundaryConnectivity_[i][j][k] << endl;
            }
        }



    }


} */

/* string MetisBoundary::SetBoundaryName(string boundaryName) {

    boundaryName_ = boundaryName;
    return boundaryName_;
} */
/*
string MetisBoundary::FindBoundaryTypeFromTagStr(string tag_str)
{

}*/