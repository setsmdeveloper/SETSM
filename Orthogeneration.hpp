//
//  Orthogeneration.hpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//

#ifndef Orthogeneration_hpp
#define Orthogeneration_hpp

#include "SubFunctions.hpp"

void orthogeneration(const TransParam _param, const ARGINFO args, char *ImageFilename, char *DEMFilename, const char *Outputpath, const int pair, const int DEM_divide, const double * const *Imageparams);
bool SetOrthoBoundary_ortho(CSize *Imagesize, double *Boundary, const double * const *RPCs, const double gridspace, const CSize DEM_size, const double minX, const double maxY, const TransParam param, const double Ortho_resolution);
bool SetDEMBoundary_ortho_photo(CSize *Imagesize, double *Boundary, const double gridspace, const CSize DEM_size, const double minX, const double maxY, const double Ortho_resolution, const EO Photo, const CAMERA_INFO m_Camera, const RM M);
uint16 *subsetImage_ortho(int check_sensor_type,FrameInfo m_frameinfo, TransParam transparam, double **RPCs, char *ImageFilename,
double *subBoundary, double *minmaxHeight, D2DPOINT *startpos, CSize* subsetsize, bool *ret);
bool GetsubareaImage_ortho(const int sensor_type, const FrameInfo m_frameinfo, const TransParam transparam, const double * const *RPCs, char *ImageFilename, CSize *Imagesize, const double *subBoundary, const double *minmaxHeight, long *cols, long *rows);
uint16 *Preprocessing_ortho(const uint8 py_level, CSize *data_size, uint16 *subimg);

#endif /* Orthogeneration_hpp */
