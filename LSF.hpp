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
CSize GetDEMsize(char *GIMP_path, char* metafilename,TransParam* param, double *grid_size, float* seeddem, double* _minX, double* _maxY);

unsigned char* GetMatchtagValue(char *GIMP_path,CSize seeddem_size);
unsigned char FloatToUnsignedChar_lsf(float val);
float UnsingedCharToFloat_lsf(float val);
void DEM_STDKenel_LSF(LSFINFO *Grid_info, double* sigma_average,double* sigma_std, float *seeddem, float *smooth_DEM, const double grid_size, const int smooth_iteration,const CSize seeddem_size, const double MPP_stereo_angle);
double LocalSurfaceFitting_DEM(LSFINFO *Grid_info, float *input, long int *numpts, double *fitted_Z, const double MPP, const int smooth_iter, const int row_size, const int col_size, const double grid, const long int X, const long int Y);
void LSFSmoothing_DEM(const char *savepath,const char* outputpath,const double MPP,const int divide);
void SetDEMBoundary_photo(EO Photo, CAMERA_INFO m_Camera, RM M, double* _boundary, double* _minmaxheight, double* _Hinterval);
#endif /* LSF_hpp */
