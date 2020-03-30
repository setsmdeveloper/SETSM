//
//  Orthogeneration.hpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//

#ifndef Orthogeneration_hpp
#define Orthogeneration_hpp

#include "SubFunctions.hpp"

void orthogeneration(TransParam _param, ARGINFO args, char *ImageFilename, char *DEMFilename, char *Outputpath,int pair,int DEM_divide,double **ImageParam);
D2DPOINT OriginalToPyramid_single_ortho(D2DPOINT InCoord, D2DPOINT Startpos, uint8 Pyramid_step);
uint16 *Preprocessing_ortho(uint8 py_level, CSize *data_size, uint16 *subimg);
uint16* CreateImagePyramid_ortho(uint16* _input, CSize _img_size, int _filter_size, double _sigma);

void SetPySizes_ortho(CSize *data_size, CSize subsetsize, int level);

uint16 *subsetImage_ortho(int check_sensor_type,FrameInfo m_frameinfo, TransParam transparam, double **RPCs, char *ImageFilename,
                          double *subBoundary, double *minmaxHeight, D2DPOINT *startpos, char *subsetImage, CSize* subsetsize, bool *ret);

bool GetImageSize_ortho(char *filename, CSize *Imagesize);
CSize Envihdr_reader_ortho(char *filename);
bool GetsubareaImage_ortho(int check_sensor_type,FrameInfo m_frameinfo, TransParam transparam, double **RPCs, char *ImageFilename, CSize *Imagesize,
                           double *subBoundary, double *minmaxHeight, int *cols, int *rows);
D2DPOINT* GetObjectToImageRPC_ortho(double **_rpc, uint8 _numofparam, double *_imageparam, uint16 _numofpts, D3DPOINT *_GP);
uint16 *Readtiff_ortho(char *filename, CSize Imagesize, int *cols, int *rows, CSize *data_size);
D2DPOINT GetObjectToImageRPC_single_ortho(double **_rpc, uint8 _numofparam, double *_imageparam, D3DPOINT _GP);


bool SetOrthoBoundary_ortho(CSize *Imagesize, double *Boundary,
                            double **RPCs, double gridspace, CSize DEM_size, double minX, double maxY, TransParam param, double Ortho_resolution);

bool SetDEMBoundary_ortho_photo(CSize *Imagesize, double *Boundary, double gridspace, CSize DEM_size, double minX, double maxY, double Ortho_resolution, EO Photo, CAMERA_INFO m_Camera, RM M);

double** OpenXMLFile_ortho(char* _filename, double* gsd_r, double* gsd_c, double* gsd);



CSize Envihdr_reader_DEM_ortho(TransParam param, char *filename, double *minX, double *maxY, double *grid_size);

#endif /* Orthogeneration_hpp */
