//
//  SubFunctions.hpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//

#ifndef SubFunctions_hpp
#define SubFunctions_hpp

#include <string.h>
#include <stdio.h>

#include <libgen.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#include <iostream>

#include "Typedefine.hpp"
#include "CoordConversion.hpp"
#include "Template.hpp"
#include "setsmgeo.hpp"
#include "grid_triangulation.hpp"
#include <sys/resource.h>
#include <sys/time.h>

#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define SUB_RATIO 4.6

char* remove_ext(const char* mystr);
char* GetFileName(char file_path[]);
char* GetFileDir(char file_path[],int *size);
char* SetOutpathName(char *_path);
int Maketmpfolders(ProInfo *info);
double getSystemMemory();

bool GetRAinfo(ProInfo *proinfo, const char* RAfile, double **Imageparams, CPairInfo &pairinfo);
bool GetRAinfoFromEcho(ProInfo *proinfo, const char* echofile, double **Imageparams, CPairInfo &pairinfo);

bool GetsubareaImage(const int sensor_type, const FrameInfo m_frameinfo, const int ti, const TransParam transparam, const double *ImageParam, const double * const *RPCs, char *ImageFilename, CSize Imagesize, const double *subBoundary, const double *minmaxHeight, long *cols, long *rows);
bool GetImageSize(char *filename, CSize *Imagesize);

bool GetsubareaImage_GeoTiff(ProInfo proinfo, char *ImageFilename, CSize Imagesize, double *subBoundary, long *cols,long *rows);

float* GetDEMValue(char *GIMP_path,CSize seeddem_size);

bool OpenProject(char* _filename, ProInfo *info, ARGINFO args);
bool OpenDMCproject(char* project_path,ProInfo *proinfo, ARGINFO args);
void SetPySizes(CSize *data_size_lr, const CSize Lsubsetsize, const int level);
void RemoveFiles(const ProInfo *proinfo,const char *save_path, char **filename, const int py_level,const bool flag);

double Correlate(const double *L, const double *R, const int N);
double Correlate(const vector<double> &L, const vector<double> &R, const int N);

void SetTranParam_fromGeoTiff(TransParam *param, char* inputfile);
void SetTransParam_param(TransParam *param, bool Hemisphere);
void SetTransParam(double minLat, double minLon, TransParam *param);

void RPCsFree(double **rpcs);
double** OpenXMLFile(ProInfo *proinfo, int ImageID, double* gsd_r, double* gsd_c, double* gsd, BandInfo* band);
double** OpenXMLFile_Pleiades(char* _filename);
double** OpenXMLFile_Planet(char* _filename);
void OpenXMLFile_orientation(char* _filename, ImageInfo *Iinfo);
void OpenXMLFile_orientation_planet(char* _filename, ImageInfo *Iinfo);

float median(int n, float* x,float min, float max);
float binmedian(int n, float *x);
double quickselect(double *arr, int n, int k);
double quickselect(vector<double> &arr, int n, int k);
bool CheckOverlap(const D2DPOINT br1_lt, const D2DPOINT br1_rb, const D2DPOINT br2_lt, const D2DPOINT br2_rb);

void Preprocessing(const ProInfo *proinfo, const char *save_path,char **Subsetfile, const uint8 py_level, const CSize * const *data_size_lr, FILE *fid);
void MakeSobelMagnitudeImage(const CSize _img_size,const uint16* _src_image, uint16* _dist_mag_image, int16* _dir);
void Orientation(const CSize imagesize,const uint16* Gmag,const int16* Gdir,const uint8 Template_size, uint8* plhs);

double InterpolatePatch(const uint16 *Image, const long int position, const CSize Imagesize, const double dx, const double dy);

struct KernelPatchArg {
    SetKernel &rkernel;

    CSize LImagesize;
    CSize RImagesize;

    const uint16 *left_image;
    const uint16 *left_mag_image;
    const uint16 *right_image;
    const uint16 *right_mag_image;
};

void SetVecKernelValue(const KernelPatchArg &kernel_patch, const int row, const int col, const D2DPOINT &pos_left, const D2DPOINT &pos_right, const int radius2, int *Count_N);

void ComputeMultiNCC(SetKernel &rsetkernel, const int Th_rho, const int *Count_N, double &count_NCC, double &sum_NCC_multi);

D2DPOINT *SetDEMGrid(const double *Boundary, const double Grid_x, const double Grid_y, CSize *Size_2D);
void SetPyramidImages(const ProInfo *proinfo, const int py_level_set, const CSize * const *data_size_lr, uint16 ***SubImages, uint16 ***SubMagImages, uint8 ***SubOriImages);
uint16 *SubsetImageFrombitsToUint16(const int image_bits, char *imagefile, long *cols, long *rows, CSize *subsize);
uint8 *SubsetImageFrombitsToUint8(const int image_bits, char *imagefile, long *cols, long *rows, CSize *subsize);

uint16* LoadPyramidImages(const char *save_path,char *subsetfile,const CSize data_size,const uint8 py_level);
uint16* LoadPyramidMagImages(const char *save_path,char *subsetfile,const CSize data_size,const uint8 py_level);
uint8* LoadPyramidOriImages(const char *save_path,char *subsetfile,const CSize data_size,const uint8 py_level);

GMA_double* GMA_double_create(uint32 size_row, uint32 size_col);
void GMA_double_destroy(GMA_double* in);
void GMA_double_inv(GMA_double *a, GMA_double *I);
void GMA_double_mul(GMA_double *a, GMA_double *b, GMA_double *out);
void GMA_double_Tran(GMA_double *a, GMA_double *out);
void GMA_double_sub(GMA_double *a, GMA_double *b, GMA_double *out);
void GMA_double_sum(GMA_double *a, GMA_double *b, GMA_double *out);
void GMA_double_printf(GMA_double *a);

FullTriangulation *TINCreate_list(D3DPOINT *ptslists, int numofpts, vector<UI3DPOINT> *trilists, double min_max[], int *count_tri, double resolution);

FullTriangulation *TINCreate_list_vector(vector<D3DPOINT> &ptslists, int numofpts, vector<UI3DPOINT> &trilists, double min_max[], int *count_tri, double resolution);

//void TINUpdate(D3DPOINT *ptslists, int numofpts, UI3DPOINT* trilists, double min_max[], int *count_tri, double resolution, FullTriangulation* oldTri, D3DPOINT* blunderlist, int numblunders);
void TINUpdate_list(D3DPOINT *ptslists, int numofpts, vector<UI3DPOINT> *trilists, double min_max[], int *count_tri, double resolution, FullTriangulation* oldTri, D3DPOINT* blunderlist, int numblunders);

void SetTinBoundary(LevelInfo &rlevelinfo, const D3DPOINT &TriP1, const D3DPOINT &TriP2, const D3DPOINT &TriP3, int *PixelMinXY, int *PixelMaxXY, double &Total_Min_Z, double &Total_Max_Z, double &temp_MinZ, double &temp_MaxZ);

bool IsTinInside(const D3DPOINT &CurGPXY, const D3DPOINT &TriP1, const D3DPOINT &TriP2, const D3DPOINT &TriP3, float &Z);

double SetNormalAngle(const D3DPOINT &pts0, const D3DPOINT &pts1, const D3DPOINT &pts2);
void SetAngle(double &angle);

//FullTriangulation *TINCreate_float(F3DPOINT *ptslists, long numofpts, UI3DPOINT* trilists, double min_max[], long *count_tri, double resolution);
//FullTriangulation *TINCreate(D3DPOINT *ptslists, int numofpts, UI3DPOINT* trilists, double min_max[], int *count_tri, double resolution);

CSize Envihdr_reader(char *filename);

CSize Envihdr_reader_seedDEM(TransParam _param, char *filename, double *minX, double *maxY, double *grid_size);
bool TFW_reader_LSFDEM(char *filename, double *minX, double *maxY, double *grid_size, int *zone, char *dir);
//bool TFW_reader_seedDEM(char *filename, double *minX, double *maxY, double *grid_size);
#endif /* SubFunctions_hpp */
