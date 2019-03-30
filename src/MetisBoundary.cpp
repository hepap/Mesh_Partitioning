#include "MetisBoundary.h"
#include <iostream>
#include <string>

MetisBoundary::MetisBoundary(int nBoundaries) 
: nBoundaries_(nBoundaries), boundaryNames_(new string[nBoundaries]), boundaryConnectivity_(nullptr), boundaryNelements_(new int[nBoundaries]), 
boundaryElementType_(new int* [nBoundaries]), boundaryElementNbrNodes_(new int* [nBoundaries]) {}


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


/* string MetisBoundary::SetBoundaryName(string boundaryName) {

    boundaryName_ = boundaryName;
    return boundaryName_;
} */
/*
string MetisBoundary::FindBoundaryTypeFromTagStr(string tag_str)
{

}*/