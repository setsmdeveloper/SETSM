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

void CollinearCalibration(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera);

void GetInitialPCfromDLT(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera);

void GetPDCs(UPARAMS param, PDCS &PD);

void EOEstimatefromInitial(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera);

void PPAEstimatefromEO(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera);

void RadialParamEstimate(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera);

void TangentialParamEstimate(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera);

void AffinityParamEstimate(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera);

void Make_b_J_K(D3DPOINT PP, D3DPOINT GP, D2DPOINT p, EO eo, CAMERA_INFO camera, CLB &B);

void CalibrationBundle(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera);
void CalibrationBundle1(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera); //all parameters 16
void CalibrationBundle2(vector<D2DPOINT> &IPs, vector<D3DPOINT> &GCPs, EO &eo, CAMERA_INFO &camera); // EO parameters 6
#endif /* CollinearEQ_hpp */
