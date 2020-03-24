#ifndef SETSM_CODE_H
#define SETSM_CODE_H

/*
 * Copyright 2017 Myoung-Jong Noh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *	   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "Typedefine.hpp"
#include "grid_triangulation.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "tiff.h"
#include "tiffio.h"
#include <list>
#include <vector>
#include <cmath>

#define	 MAXRAND	 0x7fffffff
#define	 BIGNUM		 1e37
#define	 WEEBIT		 0.000000000001
#define	 MAXDIM		 10
#define	 MAXSTR		 48
#define	 SQ(x)		 (x) * (x)
#define SWAP(a,b) temp=a;a=b;b=temp;

using std::vector;
using std::list;

typedef struct nnXY
{
    //bool flag;
	float X;
	float Y;
	float Z;
} NNXY;
void DownSample(ARGINFO &args);
//float* CreateImagePyramid_float(float* _input, CSize _img_size, int _filter_size, double _sigma);
int SETSMmainfunction(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath,double **ImageParam);
char* SetOutpathName(char *_path);

bool OpenProject(char* _filename, ProInfo *info, ARGINFO args);
int Maketmpfolders(ProInfo *info);
bool SetupParam(ProInfo *info,uint8 *NumOfIAparam, uint8 *pre_DEM_level, uint8 *DEM_level,  bool *pre_DEMtif, bool *check_tile_array);
void SetTransParam(double minLat, double minLon, bool *Hemisphere, TransParam *param);
void SetTiles(double seedDEM_gridsize, ProInfo *info, bool IsSP, bool IsRR, double *Boundary, double *Res, int tile_size, bool pre_DEMtif, uint8 *pyramid_step, uint16 *buffer_area,
			  uint8 *iter_row_start, uint8 *iter_row_end, uint8 *t_col_start, uint8 *t_col_end, double *subX, double *subY);
void SetTiles_RA(ProInfo *info, bool IsSP, bool IsRR, double *Boundary, double *Res, int tile_size, bool pre_DEMtif, uint8 *pyramid_step, uint16 *buffer_area,
				 uint8 *RA_row_start, uint8 *RA_row_end, uint8 * RA_row_iter, uint8 *t_col_start, uint8 *t_col_end, uint8 *RA_col_iter, double *subX, double *subY);
void SetPySizes(CSize *data_size_lr, const CSize Lsubsetsize, const int level);
void SetThs_ratio(int level, double *Th_roh, double *Th_roh_min, double *Th_roh_next, double *Th_roh_start, int pre_DEMtif, int IsRA, double f_demsize);

void SetThs(const ProInfo *proinfo,const int level, const int final_level_iteration, double *Th_roh, double *Th_roh_min,double *Th_roh_next, double *Th_roh_start);
D2DPOINT *SetGrids(const ProInfo *info, const int level, const int final_level_iteration, const double resolution, CSize *Size_Grid2D, const double DEM_resolution, double *py_resolution, double *grid_resolution, const double *subBoundary);
UGRID *SetGrid3PT(const ProInfo *proinfo,const TransParam param, const CSize Size_Grid2D, const double Th_roh, double *minmaxHeight, const double *subBoundary,const double py_resolution);
int	 CalTotalIteration(uint8 DEM_level,int level);

char* remove_ext(const char* mystr);

char* GetFileName(char file_path[]);
char* GetFileDir(char file_path[],int *size);

bool GetImageSize(char *filename, CSize *Imagesize);
bool GetsubareaImage(ProInfo *proinfo, int ImageID, TransParam transparam, uint8 NumofIAparam, const double * const *RPCs, const double *ImageParam, char *ImageFilename, CSize Imagesize, const double *subBoundary,const double *minmaxHeight, long int *cols,long int *rows);
uint16 *Readtiff(char *filename, CSize *Imagesize,long int *cols,long int *rows, CSize *data_size,bool check_checktiff);
bool Writetiff(char *filename, double* input, CSize Imagesize);

float *Readtiff_DEM(const char *filename, CSize *Imagesize,long int *cols,long int *rows, CSize *data_size);
unsigned char *Readtiff_BYTE(char *filename, CSize *Imagesize,long int *cols,long int *rows, CSize *data_size);
void SetSubBoundary(const double *Boundary, const double subX, const double subY, const double buffer_area, const int col, const int row, double *subBoundary);
D2DPOINT *SetDEMGrid(const double *Boundary, const double Grid_x, const double Grid_y, CSize *Size_2D);
void SetHeightWithSeedDEM(const ProInfo *proinfo,const TransParam param, UGRID *Grid, const double *Boundary, const CSize Grid_size, const double Grid_set, double *minmaxHeight);
double** OpenXMLFile(ProInfo *proinfo, int ImageID, double* gsd_r, double* gsd_c, double* gsd, BandInfo* band);
double** OpenXMLFile_Pleiades(char* _filename);
double** OpenXMLFile_Planet(char* _filename);
void OpenXMLFile_orientation(char* _filename, ImageInfo *Iinfo);

void SetDEMBoundary(double** _rpcs, double* _res,TransParam _param, bool _hemisphere, double* _boundary, double* _minmaxheight, double* _Hinterval);

bool subsetImage(ProInfo *proinfo, const TransParam transparam, const uint8 NumofIAparam, const double *const * const *RPCs, const double * const *ImageParams, const double *subBoundary, const double *minmaxHeight, D2DPOINT *Startpos, const char * const *SubsetImage, CSize *Subsetsize, FILE *fid, const bool check_checktiff);

uint16 *SetsubsetImage(ProInfo *proinfo, const int index_image, const TransParam transparam, const uint8 NumofIAparam, const double * const * const *RPCs, const double * const *ImageParams, const double *subBoundary, const double *minmaxHeight, D2DPOINT *Startpos, CSize *Subsetsize);

void SetPyramidImages(const ProInfo *proinfo, const int py_level_set, const CSize * const *data_size_lr, uint16 ***SubImages, uint16 ***SubMagImages, uint8 ***SubOriImages);

D2DPOINT* wgs2ps(TransParam _param, int _numofpts, D2DPOINT *_wgs);
D2DPOINT wgs2ps_single(TransParam _param, D2DPOINT _wgs);
D3DPOINT* wgs2ps_3D(TransParam _param, int _numofpts, D3DPOINT *_wgs);
D2DPOINT* ps2wgs(const TransParam _param,const long int _numofpts,const D2DPOINT *_ps);
D2DPOINT ps2wgs_single(TransParam _param, D2DPOINT _ps);
D2DPOINT utm2wgs_single(TransParam _param, D2DPOINT _ps);
D3DPOINT* ps2wgs_3D(TransParam _param, int _numofpts, D3DPOINT *_ps);

D2DPOINT* GetObjectToImageRPC(const double * const *_rpc, const uint8 _numofparam, const double *_imageparam, const uint16 _numofpts, D3DPOINT *_GP);
D2DPOINT GetObjectToImageRPC_single(const double * const *_rpc, const uint8 _numofparam, const double *_imageparam, const D3DPOINT _GP);
D2DPOINT GetObjectToImageRPC_single_mpp(const double * const *_rpc, const uint8 _numofparam, const double *_imageparam, const D3DPOINT _GP);
void Preprocessing(const ProInfo *proinfo, const char *save_path,char **Subsetfile, const uint8 py_level, const CSize * const *data_size_lr, FILE *fid);
uint16* LoadPyramidImages(const char *save_path,char *subsetfile,const CSize data_size,const uint8 py_level);
uint16* LoadPyramidMagImages(const char *save_path,char *subsetfile,const CSize data_size,const uint8 py_level);
uint8* LoadPyramidOriImages(const char *save_path,char *subsetfile,const CSize data_size,const uint8 py_level);

//uint16* CreateImagePyramid(uint16* _input, CSize _img_size, int _filter_size, double _sigma);
void MakeSobelMagnitudeImage(const CSize _img_size,const uint16* _src_image, uint16* _dist_mag_image, int16* _dir);
void Orientation(const CSize imagesize,const uint16* Gmag,const int16* Gdir,const uint8 Template_size, uint8* plhs);
void CalMPP_pair(double CA,double mean_product_res, double im_resolution, double *MPP_stereo_angle);
void CalMPP(const uint8 level, const CSize Size_Grid2D, const D2DPOINT* Grid_wgs,const uint8 NumofIAparam, const double * const *ImageAdjust, const double* minmaxHeight, const double * const * const *RPCs, const double CA,const double mean_product_res, const double im_resolution, double *MPP_simgle_image, double *MPP_stereo_angle);
void CalMPP_8(const CSize Size_Grid2D,const D2DPOINT* Grid_wgs,const uint8 NumofIAparam, const double * const *ImageAdjust, const double* minmaxHeight, const double * const * const *RPCs,const double CA,const double mean_product_res,const double im_resolution, double *MPP_simgle_image, double *MPP_stereo_angle);
void InitializeVoxel(const ProInfo *proinfo, VOXEL **grid_voxel,LevelInfo &plevelinfo, UGRID *GridPT3, NCCresult* nccresult,const int iteration, const double *minmaxHeight);
double GetHeightStep(int Pyramid_step, double im_resolution);

int VerticalLineLocus(VOXEL **grid_voxel,const ProInfo *proinfo, NCCresult* nccresult, LevelInfo &plevelinfo, const UGRID *GridPT3, const uint8 iteration,const double *minmaxHeight);

void SetOrthoImageCoord(const ProInfo *proinfo, LevelInfo &plevelinfo, const UGRID *GridPT3, const bool check_combined_WNCC, enum PyImageSelect check_pyimage, const double im_resolution, const double im_resolution_next, long int &sub_imagesize_w, long int &sub_imagesize_h, long int &sub_imagesize_w_next, long int &sub_imagesize_h_next, D2DPOINT **am_im_cd, D2DPOINT **am_im_cd_next);

void SetVecKernelValue(LevelInfo &plevelinfo, SetKernel &rkernel, enum PyImageSelect check_pyimage, const int row, const int col, const D2DPOINT pos_left, const D2DPOINT pos_right, const int radius2, int *Count_N);

void ComputeMultiNCC(SetKernel &rsetkernel, const int Th_rho, const int *Count_N, double &count_NCC, double &sum_NCC_multi);

void FindPeakNcc(const int Pyramid_step, const int iteration, const long int grid_index, const double temp_rho, const float iter_height, bool &check_rho, double &pre_rho, float &pre_height, int &direction, double &max_WNCC, NCCresult *nccresult);

/*
void SGM_start_pos(NCCresult *nccresult, VOXEL** grid_voxel,UGRID *GridPT3, long pt_index, bool check_ortho, double ncc_alpha, double ncc_beta, float* LHcost_pre,float **SumCost, double ortho_th, int pair_index, FILE* pfile);
void SGM_con_pos(int pts_col, int pts_row, CSize Size_Grid2D, int direction_iter, double step_height, int P_HS_step, int *u_col, int *v_row, NCCresult *nccresult, VOXEL** grid_voxel,UGRID *GridPT3, long pt_index, bool check_ortho, double ncc_alpha, double ncc_beta, double P1, double P2, float* LHcost_pre,float* LHcost_curr,float **SumCost, double ortho_th, int pair_index, FILE* pfile);
*/
void SGM_start_pos(long total_grid_size,NCCresult *nccresult, VOXEL** grid_voxel,UGRID *GridPT3, long pt_index, bool check_ortho, double ncc_alpha, double ncc_beta, float* LHcost_pre,float **SumCost, double ortho_th, int pair_index,double height_step_interval);
void SGM_con_pos(int pts_col, int pts_row, CSize Size_Grid2D, int direction_iter, double step_height, int P_HS_step, int *u_col, int *v_row, NCCresult *nccresult, VOXEL** grid_voxel,UGRID *GridPT3, long pt_index, bool check_ortho, double ncc_alpha, double ncc_beta, double P1, double P2, float* LHcost_pre,float* LHcost_curr,float **SumCost, double ortho_th, int pair_index);

void AWNCC(ProInfo *proinfo, VOXEL **grid_voxel,CSize Size_Grid2D, UGRID *GridPT3, NCCresult *nccresult, double step_height, uint8 Pyramid_step, uint8 iteration,int MaxNumberofHeightVoxel);

void rohsmoothing(double *inputroh, bool *inputcheck, int total_count, int level);

void VerticalLineLocus_seeddem(const ProInfo *proinfo,const uint16 * const *MagImages, const double DEM_resolution, const double ori_im_resolution, const double * const * const *RPCs, const CSize * const *imagesizes, const uint16 * const *Images, const uint8 Template_size, const CSize Size_Grid2D, const TransParam param, const D2DPOINT* GridPts, UGRID *GridPT3, const uint8 NumofIAparam, const double * const *ImageAdjust, const uint8 Pyramid_step, const D2DPOINT *Startpos, const double* Boundary, const double* minmaxHeight);

bool VerticalLineLocus_blunder(const ProInfo *proinfo,LevelInfo &rlevelinfo, float* nccresult, UGRID *GridPT3, uint8 iteration, bool bblunder);

int VerticalLineLocus_Ortho(ProInfo *proinfo, LevelInfo &rlevelinfo, double MPP, double *F_Height, D3DPOINT ref1_pt, D3DPOINT ref2_pt, D3DPOINT target_pt, UGRID *GridPT3, int target_index, double *F_sncc);

D2DPOINT* OriginalToPyramid(uint16 numofpts,D2DPOINT* InCoord, D2DPOINT Startpos, uint8 Pyramid_step);
D2DPOINT OriginalToPyramid_single(const D2DPOINT InCoord, const D2DPOINT Startpos, const uint8 Pyramid_step);
D2DPOINT* PyramidToOriginal(uint16 numofpts,D2DPOINT* InCoord, D2DPOINT Startpos, uint8 Pyramid_step);
long SelectMPs(const ProInfo *proinfo,LevelInfo &rlevelinfo, const NCCresult* roh_height, UGRID *GridPT3, const double Th_roh, const double Th_roh_min, const double Th_roh_start, const double Th_roh_next, const int iteration, const double MPP, const int final_level_iteration,const double MPP_stereo_angle, vector<D3DPOINT> *linkedlist);

UI3DPOINT* TINgeneration(bool last_flag, char *savepath, uint8 level, CSize Size_Grid2D, double img_resolution, double grid_resolution,
						 double min_max[],
						 double *subBoundary, int total_point_count, D3DPOINT *ptslists, int *iter_row, int *iter_col,
						 int *re_total_tri_counts);

void DecisionMPs(const ProInfo *proinfo, LevelInfo &rlevelinfo, const bool flag_blunder,const long int count_MPs_input,UGRID *GridPT3, const uint8 iteration, const double Hinterval, int *count_Results, double *minz_mp, double *maxz_mp, const double *minmaxHeight, D3DPOINT *ptslists);

void DecisionMPs_setheight(const ProInfo *proinfo, LevelInfo &rlevelinfo, const long int count_MPs_input,UGRID *GridPT3, const uint8 iteration, const double Hinterval, const double *minmaxHeight, D3DPOINT *ptslists, UI3DPOINT *trilists,int numoftri);

int SetttingFlagOfGrid(LevelInfo &rlevelinfo, UGRID *GridPT3, vector<D3DPOINT> MatchedPts_list_anchor,vector<D3DPOINT> MatchedPts_list_blunder,vector<D3DPOINT> *MatchedPts_list_mps);

int AdjustParam(ProInfo *proinfo, LevelInfo &rlevelinfo, int NumofPts, double **ImageAdjust, uint8 total_pyramid, D3DPOINT* ptslists);
bool postNCC(LevelInfo &rlevelinfo, const double Ori_diff, const D2DPOINT left_pt, const D2DPOINT right_pt, double **subA, double **TsubA, double **InverseSubA, uint8 Half_template_size, const int reference_ID, const int target_ID, double *sum_weight_X, double *sum_weight_Y, double *sum_max_roh, double **left_patch_vecs, double **right_patch_vecs);

FullTriangulation *TINCreate(D3DPOINT *ptslists, int numofpts, UI3DPOINT* trilists, double min_max[], int *count_tri, double resolution);
FullTriangulation *TINCreate_list(D3DPOINT *ptslists, int numofpts, vector<UI3DPOINT> *trilists, double min_max[], int *count_tri, double resolution);

void TINUpdate(D3DPOINT *ptslists, int numofpts, UI3DPOINT* trilists, double min_max[], int *count_tri, double resolution, FullTriangulation* oldTri, D3DPOINT* blunderlist, int numblunders);
void TINUpdate_list(D3DPOINT *ptslists, int numofpts, vector<UI3DPOINT> *trilists, double min_max[], int *count_tri, double resolution, FullTriangulation* oldTri, D3DPOINT* blunderlist, int numblunders);

bool blunder_detection_TIN(const ProInfo *proinfo, LevelInfo &rlevelinfo, const int iteration, float* ortho_ncc, bool flag_blunder, uint16 count_bl, D3DPOINT *pts, bool *detectedBlunders, long int num_points, UI3DPOINT *tris, long int num_triangles, UGRID *Gridpts, long *blunder_count,double *minz_mp, double *maxz_mp);

int Ortho_blunder(ProInfo *proinfo, LevelInfo &rlevelinfo, double MPP, D3DPOINT *pts, int numOfPts, UI3DPOINT *tris,int numOfTri, UGRID *GridPT3);

bool SetHeightRange_blunder(LevelInfo &rlevelinfo, const D3DPOINT *pts, const int numPts, UI3DPOINT *tris,const long num_triangles, UGRID *GridPT3, double *mt_minmaxheight);

UGRID* SetHeightRange(ProInfo *proinfo, LevelInfo &rlevelinfo, NCCresult *nccresult, const int numOfPts, const int num_triangles, UGRID *GridPT3, const int iteration, double *minH_grid, double *maxH_grid, D3DPOINT *pts, const UI3DPOINT *tris, const double MPP, const bool level_check_matching_rate);

UGRID* SetHeightRange_cp(ProInfo *proinfo, NCCresult *nccresult, bool pre_DEMtif, double* minmaxHeight,int numOfPts, int numOfTri, UGRID *GridPT3, bool update_flag,
                         double *minH_grid, double *maxH_grid, BL BL_param,D3DPOINT *pts, UI3DPOINT *tris,int IsRA, double MPP, char* save_path, uint8 row, uint8 col,bool check_level_end,double seedDEMsigma, bool level_check_matching_rate);
UGRID* ResizeGirdPT3(ProInfo *proinfo, CSize preSize, CSize resize_Size, double* Boundary, D2DPOINT *resize_Grid, UGRID *preGridPT3, double pre_gridsize, double* minmaxheight);
UGRID* ResizeGirdPT3_RA(const ProInfo *proinfo,const CSize preSize,const CSize resize_Size,const double* preBoundary,const double* Boundary,const D2DPOINT *resize_Grid, UGRID *preGridPT3,const double pre_gridsize,const double* minmaxheight);

void echoprint_Gridinfo(ProInfo *proinfo, NCCresult* roh_height, char *save_path,int row,int col,int level, int iteration, double update_flag, CSize *Size_Grid2D, UGRID *GridPT3, char *add_str);
void echo_print_nccresults(char *save_path,int row,int col,int level, int iteration, NCCresult *nccresult, CSize *Size_Grid2D, char *add_str);

int Matching_SETSM(ProInfo *proinfo,const uint8 pyramid_step, const uint8 Template_size, const uint16 buffer_area, const uint8 iter_row_start, const uint8 iter_row_end, const uint8 t_col_start, const uint8 t_col_end, const double subX,const double subY,const double bin_angle,const double Hinterval,const double *Image_res, double **Imageparams, const double *const*const*RPCs, const uint8 NumOfIAparam, const CSize *Imagesizes,const TransParam param, double *minmaxHeight,const double *Boundary, const double CA,const double mean_product_res, double *stereo_angle_accuracy,FILE* pMetafile);

bool check_kernel_size(ProInfo *proinfo, const CSize *Subsetsize,const int Template_size, const int pyramid_step);
bool check_image_boundary(const ProInfo *proinfo,LevelInfo &plevelinfo, const D2DPOINT pos_xy_m,const D2DPOINT pos_xy,const double minH,const double maxH,const int H_template_size);
void RemoveFiles(const ProInfo *proinfo,const char *save_path, char **filename, const int py_level,const bool flag);

short DoubleToSignedChar_result(double val);
double SignedCharToDouble_result(short val);

short DoubleToSignedChar_grid(double val);
double SignedCharToDouble_grid(short val);

short DoubleToSignedChar_voxel(double val);
double SignedCharToDouble_voxel(short val);

signed char FloatToSignedChar(float val);
float SignedCharToFloat(signed char val);

double CalMemorySize_Post(CSize DEM_size,CSize Final_DEMsize);
double CalMemorySize_Post_MT(CSize DEM_size, CSize Final_DEMsize);
double CalMemorySize_Post_LSF(CSize DEM_size, CSize Final_DEMsize);

double MergeTiles(ProInfo *info, int iter_row_start, int t_col_start, int iter_row_end,int t_col_end, int buffer,int final_iteration, float *DEM, CSize Final_DEMsize, double *FinalDEM_boundary);
double MergeTiles_Ortho(ProInfo *info, int iter_row_start, int t_col_start, int iter_row_end,int t_col_end, int buffer,int final_iteration, signed char *DEM_ortho, CSize Final_DEMsize, double *FinalDEM_boundary);

double FindNebPts_F_M_IDW(float *input, unsigned char *matching_flag, int row_size, int col_size, double grid, double minX, double minY, double maxX, double maxY, double X, double Y, int *numpts, int row_interval, int col_interval, int ndim1, char* path);

CSize DEM_final_Size(char *save_path, int row_start, int col_start,int row_end, int col_end, double grid_resolution, double *boundary);
void NNA_M(bool check_Matchtag,TransParam _param, char *save_path, char* Outputpath_name, char *iterfile, char *iterorthofile, int row_start, int col_start, int row_end, int col_end, double grid_resolution, double mt_grid_resolution, int buffer_clip, int Hemisphere,int final_iteration,int divide,CSize Final_DEMsize, float* DEM_values,float* value, unsigned char* value_pt, double *FinalDEM_boundary);
void NNA_M_MT(bool check_Matchtag,TransParam _param, char *save_path, char* Outputpath_name, char *iterfile, char *iterorthofile, int row_start, int col_start,int row_end, int col_end, double grid_resolution, double mt_grid_resolution, int buffer_clip, int Hemisphere,int final_iteration, int divide, signed char* Ortho_values,float* value, unsigned char* value_pt,CSize Final_DEMsize, double* FinalDEM_boundary);
void Envihdr_writer(TransParam _param,char *filename, int col_size, int row_size, double grid_size, double minX, double maxY, int NS_flag, int data_type);
CSize Envihdr_reader(char *filename);
CSize Envihdr_reader_seedDEM(TransParam _param, char *filename, double *minX, double *maxY, double *grid_size);
bool TFW_reader_seedDEM(char *filename, double *minX, double *maxY, double *grid_size);
bool TFW_reader_LSFDEM(char *filename, double *minX, double *maxY, double *grid_size, int *zone, char *dir);
double CalMemorySize(const ProInfo *info,LevelInfo &plevelinfo,const UGRID *GridPT3, double *minimum_memory, const uint8 iteration,const int Py_combined_level,const double *minmaxHeight);
//orthogeneration
void orthogeneration(TransParam _param, ARGINFO args, char *ImageFilename, char *DEMFilename, char *Outputpath,int pair,int DEM_divide,double **ImageParam);
D2DPOINT OriginalToPyramid_single_ortho(D2DPOINT InCoord, D2DPOINT Startpos, uint8 Pyramid_step);
uint16 *Preprocessing_ortho(uint8 py_level, CSize *data_size, uint16 *subimg);
uint16* CreateImagePyramid_ortho(uint16* _input, CSize _img_size, int _filter_size, double _sigma);

void SetPySizes_ortho(CSize *data_size, CSize subsetsize, int level);

uint16 *subsetImage_ortho(int check_sensor_type,FrameInfo m_frameinfo, TransParam transparam, double **RPCs, char *ImageFilename,
                          double *subBoundary, double *minmaxHeight, D2DPOINT *startpos, char *subsetImage, CSize* subsetsize, bool *ret);

bool GetImageSize_ortho(char *filename, CSize *Imagesize);

bool GetsubareaImage_ortho(int check_sensor_type,FrameInfo m_frameinfo, TransParam transparam, double **RPCs, char *ImageFilename, CSize *Imagesize,
                           double *subBoundary, double *minmaxHeight, int *cols, int *rows);

uint16 *Readtiff_ortho(char *filename, CSize Imagesize, int *cols, int *rows, CSize *data_size);
D2DPOINT GetObjectToImageRPC_single_ortho(double **_rpc, uint8 _numofparam, double *_imageparam, D3DPOINT _GP);
D2DPOINT* GetObjectToImageRPC_ortho(double **_rpc, uint8 _numofparam, double *_imageparam, uint16 _numofpts, D3DPOINT *_GP);

bool SetOrthoBoundary_ortho(CSize *Imagesize, double *Boundary,
                            double **RPCs, double gridspace, CSize DEM_size, double minX, double maxY, TransParam param, double Ortho_resolution);
float *LoadDEM_ortho(TransParam param, char *DEM_path, char* hdr_path);


double** OpenXMLFile_ortho(char* _filename, double* gsd_r, double* gsd_c, double* gsd);
CSize Envihdr_reader_ortho(char *filename);
CSize Envihdr_reader_DEM_ortho(TransParam param, char *filename, double *minX, double *maxY, double *grid_size);
char* remove_ext_ortho(char* mystr);


//LSF smoothing
CSize GetDEMsize(char *GIMP_path, char* metafilename,TransParam* param, double *grid_size, float* seeddem, double* _minX, double* _maxY);
float* GetDEMValue(char *GIMP_path,CSize seeddem_size);
unsigned char* GetMatchtagValue(char *GIMP_path,CSize seeddem_size);
void DEM_STDKenel_LSF(LSFINFO *Grid_info, double* sigma_average,double* sigma_std, float *seeddem, float *smooth_DEM, const double grid_size, const int smooth_iteration,const CSize seeddem_size, const double MPP_stereo_angle);
unsigned char FloatToUnsignedChar_lsf(float val);
float UnsingedCharToFloat_lsf(float val);

double LocalSurfaceFitting_DEM(LSFINFO *Grid_info, float *input, long int *numpts, double *fitted_Z, const double MPP, const int smooth_iter, const int row_size, const int col_size, const double grid, const long int X, const long int Y);
void LSFSmoothing_DEM(const char *savepath,const char* outputpath,const double MPP,const int divide);
GMA_double* GMA_double_create(uint32 size_row, uint32 size_col);
void GMA_double_destroy(GMA_double* in);
void GMA_double_inv(GMA_double *a, GMA_double *I);
void GMA_double_mul(GMA_double *a, GMA_double *b, GMA_double *out);
void GMA_double_Tran(GMA_double *a, GMA_double *out);
void GMA_double_sub(GMA_double *a, GMA_double *b, GMA_double *out);
void GMA_double_sum(GMA_double *a, GMA_double *b, GMA_double *out);
void GMA_double_printf(GMA_double *a);

RM MakeRotationMatrix(double o, double p, double k);
D2DPOINT *GetPhotoCoordinate(D3DPOINT *A, EO Photo, int _numofpts, CAMERA_INFO Camera, RM M);
D3DPOINT *GetObjectCoordinate(D2DPOINT *a, double z,EO Photo, int _numofpts, CAMERA_INFO Camera, RM M);
D2DPOINT *PhotoToImage(D2DPOINT *_photo, int _numofpts, float _CCDSize, CSize _imgsize);
D2DPOINT *ImageToPhoto(D2DPOINT *_image, int _numofpts, float _CCDSize, CSize _imgsize);

D2DPOINT GetPhotoCoordinate_single(const D3DPOINT A, const EO Photo, const CAMERA_INFO Camera, const RM M);
D3DPOINT GetObjectCoordinate_single(D2DPOINT a, double z,EO Photo, CAMERA_INFO Camera, RM M);
D2DPOINT PhotoToImage_single(const D2DPOINT _photo, const double _CCDSize, const CSize _imgsize);
D2DPOINT ImageToPhoto_single(D2DPOINT _image, float _CCDSize, CSize _imgsize);

bool OpenDMCproject(char* project_path,ProInfo *proinfo, ARGINFO args);
void SetDEMBoundary_photo(EO Photo, CAMERA_INFO m_Camera, RM M, double* _boundary, double* _minmaxheight, double* _Hinterval);
bool SetDEMBoundary_ortho_photo(CSize *Imagesize, double *Boundary, double gridspace, CSize DEM_size, double minX, double maxY, double Ortho_resolution, EO Photo, CAMERA_INFO m_Camera, RM M);
double Correlate(const double *L, const double *R, const int N);
double InterpolatePatch(const uint16 *Image, const long int position, const CSize ImageSize, const double dx, const double dy);


//Image Coregistration
double** ImageCoregistration(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath, int gcp_opt, D2DPOINT *adjust_std, bool* cal_check);
void Preprocessing_Coreg(ProInfo *proinfo, char *save_path,uint16 **Oriimage,char **Subsetfile, uint8 py_level, CSize *Subsetsize, CSize **data_size_lr);
int* CoregParam_Image(ProInfo *proinfo,uint8 Pyramid_step, uint8 total_level, double **ImageAdjust, NCCflag _flag,
                      uint8 Template_size, uint16 **Images, CSize **Imagesizes, double **Boundary, double *grid_dx, double *grid_dy,
                      int* grid_space,double** over_Boundary, char* save_path, double* avg_rho, int* iter_count, D2DPOINT *adjust_std);
bool postNCC_ortho(uint8 Pyramid_step, double Ori_diff, double Left_CR,  double Left_CC, double Right_CR, double Right_CC, double **subA,double **TsubA,double **InverseSubA, uint8 Template_size,
                   NCCflag _flag, double bin_angle, CSize leftsize, CSize rightsize, uint16* _leftimage, uint16* _rightimage, double *sum_weight_X, double *sum_weight_Y, double *sum_max_roh);
double *Readtiff_Coreg(char *filename, CSize *Imagesize, int *cols, int *rows, CSize *data_size);
double* LoadPyramidImages_double(char *save_path,char *subsetfile, CSize data_size, uint8 py_level);

//End Image Coregistration
double** DEM_ImageCoregistration(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath, int gcp_opt, D2DPOINT *adjust_std );
void DEM_ImageCoregistration_hillshade(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath, int gcp_opt);
void DEM_ImageCoregistration_GeomatricConstraint(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath, int gcp_opt);
F3DPOINT FindNormal(F3DPOINT *normal_ori, float* dem, F3DPOINT Pos,F3DPOINT Mean, F3DPOINT Scale, double* Boundary, Conformalparam X, CSize tinsize, double Gridspace, float *roh_array, float *Z, bool check_tar);
F3DPOINT SurfaceDistance_ori(F3DPOINT tar_normal_ori, float* ref_dem,F3DPOINT tar_normal, F3DPOINT tar_pts, double *tin_boundary,CSize tinsize, double Gridspace, Conformalparam param, F3DPOINT Mean, F3DPOINT Scale, float p_ref_z);
void SlopeAspect(F3DPOINT normal, F3DPOINT scale, float *slope, float *aspect);
F3DPOINT ConformalTransform(F3DPOINT input, Conformalparam param);
F3DPOINT Normalize_coord(F3DPOINT input, F3DPOINT Mean, F3DPOINT Scale);
F3DPOINT Denormalize_coord(F3DPOINT input, F3DPOINT Mean, F3DPOINT Scale);

unsigned char* CreateHillshade(float* _input, CSize _img_size, double grid_size);
F3DPOINT* SettingControls(float* DEM_ref, float* DEM_tar, double grid_size_ref, double grid_size_tar, double *boundary_ref, double *boundary_tar, double* overlapped_br, CSize img_size_ref, CSize img_size_tar, long *tin_point_num);
F3DPOINT* CreateImagePyramid_DEM(float* _input, double grid_size, double *boundary, double* overlapped_br, uint8 pyramid_level, CSize _img_size, int _filter_size, double _sigma, float* result_img, long *tin_point_num, bool check_pts);
void SetHeightRange_slope_aspect(float* ref_img, double* ref_br, CSize ref_size, float* tar_img, double* tar_br, CSize tar_size, long numOfPts, long numOfTri, F3DPOINT *pts, UI3DPOINT *tris, double *boundary, CSize input_size, double input_grid,TINinfo* tininfo, bool check_tar);
FullTriangulation *TINCreate_float(F3DPOINT *ptslists, long numofpts, UI3DPOINT* trilists, double min_max[], long *count_tri, double resolution);
D2DPOINT** CoregParam_Image_MPs(ProInfo *proinfo,uint8 Pyramid_step, uint8 total_level, double **ImageAdjust, NCCflag _flag,
                                uint8 Template_size, uint16 **Images, CSize **Imagesizes, double **Boundary, double *grid_dx, double *grid_dy,
                                int* grid_space,double** over_Boundary, char* save_path, double* avg_rho, int* iter_count, D2DPOINT *adjust_std, int*mp_iter_count);
D2DPOINT* CoregParam_Image_MPs_stereo(ProInfo *proinfo,uint8 Pyramid_step, uint8 total_level, double *ImageAdjust, NCCflag _flag,
                                      uint8 Template_size, unsigned char *Images_ref, unsigned char *Images_tar, CSize *Imagesizes_ref, CSize *Imagesizes_tar, double *Boundary_ref, double *Boundary_tar, double grid_dx_ref, double grid_dy_ref, double grid_dx_tar, double grid_dy_tar, int grid_space,double* over_Boundary, char* save_path, double* avg_rho, int* iter_count, D2DPOINT *adjust_std, int *mp_iter_count);
bool postNCC_ortho_BYTE(uint8 Pyramid_step, double Ori_diff, double Left_CR,  double Left_CC, double Right_CR, double Right_CC, double **subA,double **TsubA,double **InverseSubA, uint8 Template_size,
                        NCCflag _flag, double bin_angle, CSize leftsize, CSize rightsize, unsigned char* _leftimage, unsigned char* _rightimage, double *sum_weight_X, double *sum_weight_Y, double *sum_max_roh, F2DPOINT *peak_pos);
float SurfaceDistance(TINinfo tininfo_ref, TINinfo tininfo_tar, F3DPOINT tar_pts, F3DPOINT tar_normal, CSize tinsize, double *tin_boundary, double gridsize,Conformalparam param, int *f_iter);
float* CoeffMatrix_25D(F3DPOINT coord_center, F3DPOINT coord_scale, long pts_nums, long selected_pts, F3DPOINT* normalized_input, float* dH, Conformalparam param, F3DPOINT* tar_normal,float *weight, float* sigmaX,float *sigma0);
float* AdjustmentConformal3D(long pts_nums, long selected_pts, F3DPOINT* normalized_input, F3DPOINT* input, F3DPOINT coord_center, F3DPOINT coord_scale, float* dH, Conformalparam param, TINinfo tininfo_tar, CSize tinsize, double *tin_boundary, double gridsize, float *weight);

//SDM ortho
bool SDM_ortho(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath, double** Coreg_param);
void SetHeightWithSeedDEM_SDM(ProInfo proinfo,UGRIDSDM *Grid, double *Boundary, CSize Grid_size, double Grid_set);
int Matching_SETSM_SDM(const ProInfo proinfo,uint8 pyramid_step, uint8 Template_size, double *Image_res,double *Res, double *Limageparam, double *Rimageparam, CSize Limagesize,CSize Rimagesize, double *Boundary, ImageGSD gsd_image1, ImageGSD gsd_image2,double* LBoundary,double* RBoundary, int *matching_number);
bool subsetImage_SDM(ProInfo proinfo,char *LImageFilename, char *RImageFilename, double *subBoundary, D2DPOINT *Lstartpos, D2DPOINT *Rstartpos, char *LsubsetImage, char *RsubsetImage, CSize* Lsubsetsize, CSize* Rsubsetsize, FILE *fid,bool check_checktiff);
bool GetsubareaImage_SDM(ProInfo proinfo, char *ImageFilename, CSize *Imagesize, double *subBoundary, long int *cols,long int *rows);
D2DPOINT* GetObjectToImage(uint16 _numofpts, D2DPOINT *_GP, double *boundary, double imageres);
void Preprocessing_SDM(ProInfo proinfo, char *save_path,char *Lsubsetfile, char *Rsubsetfile, uint8 py_level, CSize *Lsubsetsize, CSize *Rsubsetsize, CSize *data_size_l, CSize *data_size_r, FILE *fid);
void SetThs_SDM(int level, int final_level_iteration, double *Th_roh, double *Th_roh_min, double *Th_roh_next, double *Th_roh_start, uint8 pyramid_step);
D2DPOINT *SetGrids_SDM(ProInfo proinfo, int prc_level,int level, int start_py, int final_level_iteration, double resolution, CSize *Size_Grid2D, double DEM_resolution, double *py_resolution, double *grid_resolution, double *subBoundary);
D2DPOINT *SetDEMGrid_SDM(double *Boundary, double Grid_x, double Grid_y, CSize *Size_2D);
UGRIDSDM *SetGrid3PT_SDM(ProInfo proinfo,bool flag_start, CSize Size_Grid2D, double Th_roh, double *subBoundary, double py_resolution);
UGRIDSDM* ResizeGirdPT3_SDM(CSize preSize, CSize resize_Size, double* Boundary, D2DPOINT *resize_Grid, UGRIDSDM *preGridPT3, double pre_gridsize);
bool VerticalLineLocus_SDM(ProInfo proinfo, NCCresultSDM* nccresult, uint16 *MagImages_L,uint16 *MagImages_R,double DEM_resolution, double im_resolution, CSize LImagesize_ori, CSize LImagesize, uint16* LeftImage, CSize RImagesize_ori, CSize RImagesize, uint16* RightImage, uint8 Template_size, CSize Size_Grid2D, D2DPOINT* GridPts, UGRIDSDM *GridPT3, uint8 Pyramid_step,uint8 start_py, D2DPOINT Lstartpos, D2DPOINT Rstartpos, uint8 iteration, double* Boundary, ImageGSD gsd_image1, ImageGSD gsd_image2, double* Coreg_param,uint8* SubOriImages_L,uint8* SubOriImages_R);
int SelectMPs_SDM(ProInfo proinfo, NCCresultSDM* roh_height, CSize Size_Grid2D, D2DPOINT *GridPts_XY, UGRIDSDM *GridPT3, uint8 Pyramid_step,uint8 prc_level, uint8 iteration, char *filename_mps,uint8 tile_row,uint8 tile_col,double *Boundary);
void echo_print_nccresults_SDM(char *save_path,int row,int col,int level, int iteration, NCCresultSDM *nccresult, CSize *Size_Grid2D, char *add_str);
void echoprint_Gridinfo_SDM(uint8 prc_level, ProInfo proinfo, double* LBoundary,double* RBoundary, CSize LImagesize, uint16* LeftImage, CSize RImagesize, uint16* RightImage,double* boundary,char *save_path,int row,int col,int level, int iteration, CSize *Size_Grid2D, UGRIDSDM *GridPT3, const char *add_str);
void echoprint_adjustXYZ(uint8 prc_level, double* LBoundary,double* RBoundary, CSize LImagesize, uint16* LeftImage, CSize RImagesize, uint16* RightImage,double* boundary,ProInfo proinfo, char *save_path,int row,int col,int level, int iteration, CSize *Size_Grid2D, UGRIDSDM *GridPT3,const char *add_str, double gridsize, int d_date);
UGRIDSDM* SetHeightRange_SDM(UGRIDSDM *GridPT3, BL BL_param);
void SetHeightRange_shift(int numOfPts, int numOfTri, UGRIDSDM *GridPT3,BL BL_param,D3DPOINT *pts, UI3DPOINT *tris,int b_dir);
bool Update_ortho_NCC(ProInfo proinfo, uint16 *MagImages_L,uint16 *MagImages_R,double DEM_resolution, double im_resolution, CSize LImagesize_ori, CSize LImagesize, uint16* LeftImage, CSize RImagesize_ori, CSize RImagesize, uint16* RightImage, uint8 Template_size, CSize Size_Grid2D, D2DPOINT* GridPts, UGRIDSDM *GridPT3, uint8 Pyramid_step, D2DPOINT Lstartpos, D2DPOINT Rstartpos, uint8 iteration, double* Boundary, ImageGSD gsd_image1, ImageGSD gsd_image2, double* Coreg_param,uint8* SubOriImages_L,uint8* SubOriImages_R);
bool average_filter_colrowshift(CSize Size_Grid2D, UGRIDSDM *GridPT3,uint8 Pyramid_step);
void shift_filtering(ProInfo proinfo,UGRIDSDM *GridPT3, BL BL_param, double DEM_resolution);
void RemoveFiles_SDM(char *save_path, char *lfilename, char *rfilename, int py_level, bool flag);
double MergeTiles_SDM(ProInfo info,int iter_row_end,int t_col_end, int buffer,int final_iteration,TransParam _param,int Hemisphere,uint8 pyramid_step);
D2DPOINT GetObjectToImage_single_SDM(uint16 _numofpts, D2DPOINT _GP, double *boundary, double imageres);
void SetTransParam_param(TransParam *param, bool Hemisphere);
void SetTranParam_fromOrtho(TransParam *param,char* inputfile,ARGINFO args,bool *Hemisphere);
float median(int n, float* x,float min, float max);
float binmedian(int n, float *x);
float quickselect(float *arr, int n, int k);
bool CheckOverlap(D2DPOINT br1_lt, D2DPOINT br1_rb, D2DPOINT br2_lt, D2DPOINT br2_rb);

inline double SQRT(D2DPOINT a);
inline double SQRT(D2DPOINT a, D2DPOINT b);
inline double SQRT(D3DPOINT a, int dimension = 3);
inline double SQRT(D3DPOINT a, D3DPOINT b, int dimension = 3);

template <typename T>
T* CreateImagePyramid(T* _input, CSize _img_size, int _filter_size, double _sigma);
template <typename T>
T BilinearResampling(T* input, const CSize img_size, D2DPOINT query_pt);
#endif // SETSM_CODE_H
