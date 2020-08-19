//
//  LSF.hpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//

#ifndef LSF_hpp
#define LSF_hpp

#include "SubFunctions.hpp"

//LSF smoothing
void LSFSmoothing_DEM(const char *savepath,const char* outputpath,const double MPP,const int divide);
CSize GetDEMsize(char *GIMP_path, char* metafilename,TransParam* param, double *grid_size, double* _minX, double* _maxY);
void DEM_STDKenel_LSF(LSFINFO *Grid_info, double* sigma_average,double* sigma_std, float *seeddem, float *smooth_DEM, const double grid_size, const int smooth_iteration,const CSize seeddem_size, const double MPP_stereo_angle);
double LocalSurfaceFitting_DEM(LSFINFO *Grid_info, float *input, long &numpts, double *fitted_Z, const double MPP, const int smooth_iter, const long row_size, const long col_size, const double grid, const long X, const long Y);
#endif /* LSF_hpp */
