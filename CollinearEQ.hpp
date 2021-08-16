//
//  CollinearEQ.hpp
//  
//
//  Created by MJ_Noh on 8/6/21.
//

#ifndef CollinearEQ_hpp
#define CollinearEQ_hpp

#include <cmath>
#include <stdlib.h>

#include "Typedefine.hpp"
#include "SubFunctions.hpp"

void GetInitialPCfromDLT(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera);

void GetPDCs(UPARAMS param, PDCS &PD);

void CalibrationBundle(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera);

#endif /* CollinearEQ_hpp */
