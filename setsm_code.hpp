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

#define	 MAXRAND	 0x7fffffff
#define	 BIGNUM		 1e37
#define	 WEEBIT		 0.000000000001
#define	 MAXDIM		 10
#define	 MAXSTR		 48
#define	 SQ(x)		 (x) * (x)
#define SWAP(a,b) temp=a;a=b;b=temp;

typedef struct nnXY
{
    //bool flag;
	float X;
	float Y;
	float Z;
} NNXY;

int SETSMmainfunction(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath,double **ImageParam);
char* SetOutpathName(char *_path);

bool OpenProject(char* _filename, ProInfo *info, ARGINFO args);
int Maketmpfolders(ProInfo *info);
bool SetupParam(ProInfo *info,uint8 *NumOfIAparam, uint8 *pre_DEM_level, uint8 *DEM_level,  bool *pre_DEMtif, bool *check_tile_array);
void SetTransParam(double minLat, double minLon, bool *Hemisphere, TransParam *param);
void SetTiles(ProInfo *info, bool IsSP, bool IsRR, double *Boundary, double *Res, int tile_size, bool pre_DEMtif, uint8 *pyramid_step, uint16 *buffer_area,
			  uint8 *iter_row_start, uint8 *iter_row_end, uint8 *t_col_start, uint8 *t_col_end, double *subX, double *subY);
void SetTiles_RA(ProInfo *info, bool IsSP, bool IsRR, double *Boundary, double *Res, int tile_size, bool pre_DEMtif, uint8 *pyramid_step, uint16 *buffer_area,
				 uint8 *RA_row_start, uint8 *RA_row_end, uint8 * RA_row_iter, uint8 *t_col_start, uint8 *t_col_end, uint8 *RA_col_iter, double *subX, double *subY);
void SetPySizes(CSize *data_size_lr, CSize Lsubsetsize, int level);
void SetThs_ratio(int level, double *Th_roh, double *Th_roh_min, double *Th_roh_next, double *Th_roh_start, int pre_DEMtif, int IsRA, double f_demsize);
void SetThs(ProInfo *proinfo,int level, int final_level_iteration, double *Th_roh, double *Th_roh_min, double *Th_roh_next, double *Th_roh_start);
D2DPOINT *SetGrids(ProInfo *info, bool *dem_update_flag, bool flag_start, int level, int final_level_iteration, double resolution, CSize *Size_Grid2D, bool pre_DEMtif, char *priori_DEM_tif, double DEM_resolution, double *minmaxHeight,
				   double *py_resolution, double *grid_resolution, double *subBoundary);
UGRID *SetGrid3PT(ProInfo *proinfo,TransParam param, bool dem_update_flag, bool flag_start, CSize Size_Grid2D, double Th_roh, int level, double *minmaxHeight,double *subBoundary,double py_resolution,char* metafilename);
int	 CalTotalIteration(uint8 DEM_level,int level);

char* remove_ext(char* mystr);

char* GetFileName(char file_path[]);
char* GetFileDir(char file_path[],int *size);

bool GetImageSize(char *filename, CSize *Imagesize);
bool GetsubareaImage(ProInfo *proinfo, int ImageID, TransParam transparam, uint8 NumofIAparam, double **RPCs, double *ImageParam, char *ImageFilename, CSize Imagesize,
					 double *subBoundary, double *minmaxHeight, int *cols, int *rows);
uint16 *Readtiff(char *filename, CSize *Imagesize, int *cols, int *rows, CSize *data_size,bool check_checktiff);
bool Writetiff(char *filename, double* input, CSize Imagesize);

float *Readtiff_DEM(char *filename, CSize *Imagesize, int *cols, int *rows, CSize *data_size);
unsigned char *Readtiff_BYTE(char *filename, CSize *Imagesize, int *cols, int *rows, CSize *data_size);
void SetSubBoundary(double *Boundary, double subX, double subY, double buffer_area, int col, int row, double *subBoundary);
D2DPOINT *SetDEMGrid(double *Boundary, double Grid_x, double Grid_y, CSize *Size_2D);
void SetHeightWithSeedDEM(ProInfo *proinfo,TransParam param, UGRID *Grid, double *Boundary, CSize Grid_size, double Grid_set, double *minmaxHeight);

double** OpenXMLFile(ProInfo *proinfo, int ImageID, double* gsd_r, double* gsd_c, double* gsd, BandInfo* band);
double** OpenXMLFile_Pleiades(char* _filename);
double** OpenXMLFile_Planet(char* _filename);
void OpenXMLFile_orientation(char* _filename, ImageInfo *Iinfo);

void SetDEMBoundary(double** _rpcs, double* _res,TransParam _param, bool _hemisphere, double* _boundary, double* _minmaxheight, double* _Hinterval);

bool subsetImage(ProInfo *proinfo, TransParam transparam, uint8 NumofIAparam, double ***RPCs, double **ImageParams,
				 double *subBoundary, double *minmaxHeight, D2DPOINT *Startpos, char **SubsetImage, CSize *Subsetsize, FILE *fid,bool check_checktiff);
D2DPOINT* wgs2ps(TransParam _param, int _numofpts, D2DPOINT *_wgs);
D2DPOINT wgs2ps_single(TransParam _param, D2DPOINT _wgs);
D3DPOINT* wgs2ps_3D(TransParam _param, int _numofpts, D3DPOINT *_wgs);
D2DPOINT* ps2wgs(TransParam _param, int _numofpts, D2DPOINT *_ps);
D2DPOINT ps2wgs_single(TransParam _param, D2DPOINT _ps);
D2DPOINT utm2wgs_single(TransParam _param, D2DPOINT _ps);
D3DPOINT* ps2wgs_3D(TransParam _param, int _numofpts, D3DPOINT *_ps);

D2DPOINT* GetObjectToImageRPC(double **_rpc, uint8 _numofparam, double *_imageparam, uint16 _numofpts, D3DPOINT *_GP);
D2DPOINT GetObjectToImageRPC_single(double **_rpc, uint8 _numofparam, double *_imageparam, D3DPOINT _GP);
D2DPOINT GetObjectToImageRPC_single_mpp(double **_rpc, uint8 _numofparam, double *_imageparam, D3DPOINT _GP);
void Preprocessing(ProInfo *proinfo, char *save_path,char **Subsetfile, uint8 py_level, CSize *Subsetsize, CSize **data_size_lr, FILE *fid);
uint16* LoadPyramidImages(char *save_path,char *subsetfile, CSize data_size, uint8 py_level);
uint16* LoadPyramidMagImages(char *save_path,char *subsetfile, CSize data_size, uint8 py_level, double *val, double *avg);
uint8* LoadPyramidOriImages(char *save_path,char *subsetfile, CSize data_size, uint8 py_level);
uint16* CreateImagePyramid(uint16* _input, CSize _img_size, int _filter_size, double _sigma);
void MakeSobelMagnitudeImage(CSize _img_size, uint16* _src_image, uint16* _dist_mag_image, /*int16* _gx, int16* _gy,*/ int16* _dir);
void Orientation(CSize imagesize, uint16* Gmag, int16* Gdir, uint8 Template_size, uint8* plhs);
void CalMPP_pair(double CA,double mean_product_res, double im_resolution, double *MPP_stereo_angle);
void CalMPP(uint8 level, CSize Size_Grid2D, TransParam param, D2DPOINT* Grid_wgs,uint8 NumofIAparam, double **ImageAdjust, double* minmaxHeight, double ***RPCs, double CA,double mean_product_res, double im_resolution, double *MPP_simgle_image, double *MPP_stereo_angle);
void CalMPP_8(uint8 level, CSize Size_Grid2D, TransParam param, D2DPOINT* Grid_wgs,uint8 NumofIAparam, double **ImageAdjust, double* minmaxHeight, double ***RPCs,double CA,double mean_product_res, double im_resolution, double *MPP_simgle_image, double *MPP_stereo_angle);
void InitializeVoxel(VOXEL **grid_voxel,CSize Size_Grid2D, double height_step, UGRID *GridPT3, NCCresult* nccresult, int iteration, uint8 pyramid_step, double DEM_resolution,ProInfo *proinfo,double ***RPCs,double **ImageAdjust,D2DPOINT *Startpos,D2DPOINT* GridPts, D2DPOINT* Grid_wgs,CSize *Imagesizes_ori, CSize **Imagesizes,double* minmaxHeight);

double GetHeightStep(int Pyramid_step, double im_resolution);
int VerticalLineLocus(VOXEL **grid_voxel, ProInfo *proinfo, NCCresult* nccresult, uint16 **MagImages,double DEM_resolution, double im_resolution, double ***RPCs, CSize *Imagesizes_ori, CSize **Imagesizes, uint16 **Images, uint8 Template_size,
					   CSize Size_Grid2D, TransParam param, D2DPOINT* GridPts, D2DPOINT* Grid_wgs, UGRID *GridPT3, NCCflag flag,
					   uint8 NumofIAparam, double **ImageAdjust, double* minmaxHeight, uint8 Pyramid_step, D2DPOINT *Startpos, uint8 iteration, uint8 **ori_images,
					   double bin_angle, uint8 NumOfCompute, uint8 peak_level, FILE* fid, bool IsPar, bool Hemisphere, uint8 tile_row, uint8 tile_col, double* Boundary,
					   char* tmpdir, double mag_avg,double mag_var,D2DPOINT *Startpos_next,uint16 **SubImages_next,uint8 **SubOriImages_next,uint16 **SubMagImages_next,int Py_combined_level, bool check_matching_rate);
/*
void SGM_start_pos(NCCresult *nccresult, VOXEL** grid_voxel,UGRID *GridPT3, long pt_index, bool check_ortho, double ncc_alpha, double ncc_beta, float* LHcost_pre,float **SumCost, double ortho_th, int pair_index, FILE* pfile);
void SGM_con_pos(int pts_col, int pts_row, CSize Size_Grid2D, int direction_iter, double step_height, int P_HS_step, int *u_col, int *v_row, NCCresult *nccresult, VOXEL** grid_voxel,UGRID *GridPT3, long pt_index, bool check_ortho, double ncc_alpha, double ncc_beta, double P1, double P2, float* LHcost_pre,float* LHcost_curr,float **SumCost, double ortho_th, int pair_index, FILE* pfile);
*/
void SGM_start_pos(long total_grid_size,NCCresult *nccresult, VOXEL** grid_voxel,UGRID *GridPT3, long pt_index, bool check_ortho, double ncc_alpha, double ncc_beta, float* LHcost_pre,float **SumCost, double ortho_th, int pair_index);
void SGM_con_pos(int pts_col, int pts_row, CSize Size_Grid2D, int direction_iter, double step_height, int P_HS_step, int *u_col, int *v_row, NCCresult *nccresult, VOXEL** grid_voxel,UGRID *GridPT3, long pt_index, bool check_ortho, double ncc_alpha, double ncc_beta, double P1, double P2, float* LHcost_pre,float* LHcost_curr,float **SumCost, double ortho_th, int pair_index);

void AWNCC(ProInfo *proinfo, VOXEL **grid_voxel,CSize Size_Grid2D, UGRID *GridPT3, NCCresult *nccresult, double step_height, uint8 Pyramid_step, uint8 iteration,int MaxNumberofHeightVoxel);

void rohsmoothing(double *inputroh, bool *inputcheck, int total_count, int level);

double VerticalLineLocus_seeddem(ProInfo *proinfo,uint16 **MagImages, double DEM_resolution, double im_resolution, double ***RPCs,
								CSize *Imagesizes_ori, CSize **Imagesizes, uint16 **Images, uint8 Template_size,
								CSize Size_Grid2D, TransParam param, D2DPOINT* GridPts, D2DPOINT *Grid_wgs, UGRID *GridPT3,
								uint8 NumofIAparam, double **ImageAdjust, uint8 Pyramid_step, D2DPOINT *Startpos,
								char* save_filepath, uint8 tile_row, uint8 tile_col, uint8 iteration,uint8 bl_count,double* Boundary, double* minmaxHeight, double seedDEMsigma);

bool VerticalLineLocus_blunder(ProInfo *proinfo,double* nccresult, double* INCC, uint16 **MagImages, double DEM_resolution, double im_resolution, double ***RPCs,
                               CSize *Imagesizes_ori, CSize **Imagesizes, uint16 **Images, uint8 Template_size,
                               CSize Size_Grid2D, TransParam param, D2DPOINT* GridPts, D2DPOINT *Grid_wgs, UGRID *GridPT3,
                               uint8 NumofIAparam, double **ImageAdjust, uint8 Pyramid_step, D2DPOINT *Startpos,
                               char* save_filepath, uint8 tile_row, uint8 tile_col, uint8 iteration,uint8 bl_count,double* Boundary,uint8 **ori_images, int blunder_selected_level, bool bblunder);
int VerticalLineLocus_Ortho(ProInfo *proinfo, double *F_Height,D3DPOINT ref1_pt, D3DPOINT ref2_pt, D3DPOINT target_pt,
                             uint16 **MagImages, uint16 **Images,
                             double DEM_resolution, double im_resolution, double ***RPCs,
                             CSize **Imagesizes,  CSize Size_Grid2D, TransParam param, uint8 NumofIAparam,
                             double **ImageAdjust, double* minmaxHeight, uint8 Pyramid_step, double meters_per_pixel,
                             D2DPOINT *Startpos, uint8 iteration,  UGRID *GridPT3, int target_index, int ref1_index, int ref2_index,
                             double* boundary,double *F_sncc);
D2DPOINT* OriginalToPyramid(uint16 numofpts,D2DPOINT* InCoord, D2DPOINT Startpos, uint8 Pyramid_step);
D2DPOINT OriginalToPyramid_single(D2DPOINT InCoord, D2DPOINT Startpos, uint8 Pyramid_step);
D2DPOINT* PyramidToOriginal(uint16 numofpts,D2DPOINT* InCoord, D2DPOINT Startpos, uint8 Pyramid_step);
int SelectMPs(ProInfo *proinfo,NCCresult* roh_height, CSize Size_Grid2D, D2DPOINT *GridPts_XY, UGRID *GridPT3,
			  double Th_roh, double Th_roh_min, double Th_roh_start, double Th_roh_next, uint8 Pyramid_step, uint8 total_pyramid,
			  uint8 iteration, uint8 peak_level, char *filename_mps, int pre_DEMtif, int IsRA, double MPP, double DEM_resolution, double im_resolution, int final_level_iteration,double MPP_stereo_angle);

UI3DPOINT* TINgeneration(bool last_flag, char *savepath, uint8 level, CSize Size_Grid2D, double img_resolution, double grid_resolution,
						 double min_max[],
						 double *subBoundary, int total_point_count, D3DPOINT *ptslists, int *iter_row, int *iter_col,
						 int *re_total_tri_counts);

int DecisionMPs(ProInfo *proinfo,bool flag_blunder,int count_MPs, double* Boundary, UGRID *GridPT3, uint8 Pyramid_step, double grid_resolution,
                uint8 iteration, CSize Size_Grid2D, char *filename_mps_pre, char *filename_mps, double Hinterval,
                bool *p_flag, double *pre_3sigma, double *pre_mean, int *count_Results, double *minz_mp, double *maxz_mp, double *minmaxHeight,
                uint16 **MagImages,double DEM_resolution, double im_resolution, double ***RPCs,
                CSize *Imagesizes_ori, CSize **Imagesizes, uint16 **Images, uint8 Template_size,
                TransParam param, D2DPOINT* Grid_wgs,D2DPOINT* GridPts,
                uint8 NumofIAparam, double **ImageAdjust, D2DPOINT *Startpos,
                uint8 tile_row, uint8 tile_col, uint8 **ori_images, int blunder_selected_level);

int DecisionMPs_setheight(ProInfo *proinfo,bool flag_blunder, int count_MPs_input, double* Boundary, UGRID *GridPT3, uint8 Pyramid_step, double grid_resolution,
						  uint8 iteration, CSize Size_Grid2D, char *filename_mps_pre, char *filename_tri, double Hinterval,
						  bool *p_flag, double *pre_3sigma, double *pre_mean, int *count_Results, double *minz_mp, double *maxz_mp, double *minmaxHeight,
						  uint16 **MagImages,double DEM_resolution, double im_resolution, double ***RPCs,
						  CSize *Imagesizes_ori, CSize **Imagesizes, uint16 **Images, uint8 Template_size,
						  TransParam param, D2DPOINT* Grid_wgs,D2DPOINT* GridPts,
						  uint8 NumofIAparam, double **ImageAdjust, D2DPOINT *Startpos,
						  char* save_filepath, uint8 tile_row, uint8 tile_col,D3DPOINT *ptslists, UI3DPOINT *trilists,int numoftri,uint8 **ori_images, int blunder_selected_level);
int SetttingFlagOfGrid(double *subBoundary,UGRID *GridPT3, uint8 Pyramid_step,double grid_resolution,uint8 iteration,
					   CSize Size_Grid2D,char *filename_mps_anchor,char *filename_mps_aft,int count_results_anchor,int count_results_blunder, char *filename_mps);

int AdjustParam(ProInfo *proinfo,uint8 Pyramid_step, int NumofPts, char * file_pts, D2DPOINT *Startpos, double ***RPCs, double **ImageAdjust, NCCflag _flag,
				uint8 Template_size, uint16 **Images, CSize **Imagesizes, uint8 **ori_images, TransParam param,
				double bin_angle, uint8 total_pyramid, bool Hemisphere, char* save_filepath, char* tmpdir);
bool postNCC(uint8 Pyramid_step, double Ori_diff, double Left_CR,  double Left_CC, double Right_CR, double Right_CC, double **subA,double **TsubA,double **InverseSubA, uint8 Half_template_size, 
			 NCCflag _flag, double bin_angle, CSize leftsize, CSize rightsize, uint16* _leftimage, uint16* _rightimage, double *sum_weight_X, double *sum_weight_Y, double *sum_max_roh,double **left_patch_vecs,double **right_patch_vecs);

FullTriangulation *TINCreate(D3DPOINT *ptslists, int numofpts, UI3DPOINT* trilists, double min_max[], int *count_tri, double resolution);
void TINUpdate(D3DPOINT *ptslists, int numofpts, UI3DPOINT* trilists, double min_max[], int *count_tri, double resolution, FullTriangulation* oldTri, D3DPOINT* blunderlist, int numblunders);
bool blunder_detection_TIN(int pre_DEMtif, double* ortho_ncc, double* INCC,bool flag_blunder,uint16 count_bl,double* blunder_dh,char *file_pts,
						   D3DPOINT *ptslists, bool *detectedBlunders, int numOfPts, UI3DPOINT *trilists,int numOfTri, UGRID *Gridpts, BL BL_param, 
						   uint32 *blunder_count,double *minz_mp, double *maxz_mp, double *minmaxHeight, int IsRA,double seedDEMsigma);

int Ortho_blunder(ProInfo *proinfo, D3DPOINT *pts, int numOfPts, UI3DPOINT *tris,int numOfTri, bool update_flag,double *minH_grid, double *maxH_grid, BL BL_param,
				  uint16 **MagImages, uint16 **Images,
				  double DEM_resolution, double im_resolution, double ***RPCs,
				  CSize **Imagesizes, CSize Size_Grid2D, TransParam param, uint8 NumofIAparam,
				  double **ImageAdjust, double* minmaxHeight, uint8 Pyramid_step, double meters_per_pixel,
				  D2DPOINT *Startpos, uint8 iteration,	UGRID *GridPT3, char *filename_mps);

bool SetHeightRange_blunder(double* minmaxHeight,D3DPOINT *pts, int numOfPts, UI3DPOINT *tris,int numOfTri, UGRID *GridPT3, BL BL_param, double *mt_minmaxheight,bool blunder_update);
UGRID* SetHeightRange(ProInfo *proinfo, NCCresult *nccresult, bool pre_DEMtif, double* minmaxHeight,int numOfPts, int numOfTri, UGRID *GridPT3, bool update_flag,
					  double *minH_grid, double *maxH_grid, BL BL_param,D3DPOINT *pts, UI3DPOINT *tris,int IsRA, double MPP, char* save_path, uint8 tile_row, uint8 tile_col,bool check_level_end, double seedDEMsigma, bool level_check_matching_rate);
UGRID* SetHeightRange_cp(ProInfo *proinfo, NCCresult *nccresult, bool pre_DEMtif, double* minmaxHeight,int numOfPts, int numOfTri, UGRID *GridPT3, bool update_flag,
                         double *minH_grid, double *maxH_grid, BL BL_param,D3DPOINT *pts, UI3DPOINT *tris,int IsRA, double MPP, char* save_path, uint8 row, uint8 col,bool check_level_end,double seedDEMsigma, bool level_check_matching_rate);
UGRID* ResizeGirdPT3(ProInfo *proinfo, CSize preSize, CSize resize_Size, double* Boundary, D2DPOINT *resize_Grid, UGRID *preGridPT3, double pre_gridsize, double* minmaxheight);
UGRID* ResizeGirdPT3_RA(ProInfo *proinfo, CSize preSize, CSize resize_Size, double* preBoundary,double* Boundary, D2DPOINT *resize_Grid, UGRID *preGridPT3, double pre_gridsize, double* minmaxheight);

void echoprint_Gridinfo(ProInfo *proinfo, NCCresult* roh_height, char *save_path,int row,int col,int level, int iteration, double update_flag, CSize *Size_Grid2D, UGRID *GridPT3, char *add_str);
void echo_print_nccresults(char *save_path,int row,int col,int level, int iteration, NCCresult *nccresult, CSize *Size_Grid2D, char *add_str);

int Matching_SETSM(ProInfo *proinfo,uint8 pyramid_step, uint8 Template_size, uint16 buffer_area,uint8 iter_row_start, uint8 iter_row_end,uint8 t_col_start,uint8 t_col_end,
				   double subX,double subY,double bin_angle,double Hinterval,double *Image_res,double *Res, double *Limageparam, double **Rimageparam,
				   double ***RPCs, uint8 pre_DEM_level, uint8 DEM_level,	uint8 NumOfIAparam, bool check_tile_array,bool Hemisphere,bool* tile_array,
				   CSize *Imagesizes,TransParam param,int total_count,double *ori_minmaxHeight,double *Boundary,int row_iter, int col_iter, double CA,double mean_product_res, double *stereo_angle_accuracy,FILE* pMetafile);
bool check_kernel_size(ProInfo *proinfo, CSize *Subsetsize, int Template_size, int pyramid_step);
bool check_image_boundary(ProInfo *proinfo, double ***rpc, uint8 numofparam, double **imageparam, D2DPOINT *Startpos,
						  D2DPOINT pos_xy_m,D2DPOINT pos_xy, double minH, double maxH, CSize *Imagesizes_ori, CSize **sizes, int H_template_size, int pyramid_step);
void RemoveFiles(ProInfo *proinfo,char *save_path, char **filename, int py_level, bool flag);

double CalMemorySize_Post(CSize DEM_size,CSize Final_DEMsize);
double CalMemorySize_Post_MT(CSize DEM_size, CSize Final_DEMsize);
double CalMemorySize_Post_LSF(CSize DEM_size, CSize Final_DEMsize);

double MergeTiles(ProInfo *info, int iter_row_start, int t_col_start, int iter_row_end,int t_col_end, int buffer,int final_iteration, float *DEM, CSize Final_DEMsize, double *FinalDEM_boundary);
double MergeTiles_Ortho(ProInfo *info, int iter_row_start, int t_col_start, int iter_row_end,int t_col_end, int buffer,int final_iteration, float *DEM_ortho, CSize Final_DEMsize, double *FinalDEM_boundary);

double FindNebPts_F_M_IDW(float *input, unsigned char *matching_flag, int row_size, int col_size, double grid, double minX, double minY, double maxX, double maxY, double X, double Y, int *numpts, int row_interval, int col_interval, int ndim1, char* path);

CSize DEM_final_Size(char *save_path, int row_start, int col_start,int row_end, int col_end, double grid_resolution, double *boundary);
void NNA_M(bool check_Matchtag,TransParam _param, char *save_path, char* Outputpath_name, char *iterfile, char *iterorthofile, int row_start, int col_start, int row_end, int col_end, double grid_resolution, double mt_grid_resolution, int buffer_clip, int Hemisphere,int final_iteration,int divide,CSize Final_DEMsize, float* DEM_values,float* value, unsigned char* value_pt, double *FinalDEM_boundary);
void NNA_M_MT(bool check_Matchtag,TransParam _param, char *save_path, char* Outputpath_name, char *iterfile, char *iterorthofile, int row_start, int col_start,int row_end, int col_end, double grid_resolution, double mt_grid_resolution, int buffer_clip, int Hemisphere,int final_iteration, int divide, float* Ortho_values,float* value, unsigned char* value_pt,CSize Final_DEMsize, double* FinalDEM_boundary);
void Envihdr_writer(TransParam _param,char *filename, int col_size, int row_size, double grid_size, double minX, double maxY, int NS_flag, int data_type);
CSize Envihdr_reader(char *filename);
CSize Envihdr_reader_seedDEM(TransParam _param, char *filename, double *minX, double *maxY, double *grid_size);
bool TFW_reader_seedDEM(char *filename, double *minX, double *maxY, double *grid_size);
bool TFW_reader_LSFDEM(char *filename, double *minX, double *maxY, double *grid_size, int *zone, char *dir);
double CalMemorySize(ProInfo *info, CSize Size_Grid2D,CSize** data_size,UGRID *GridPT3, int level, double height_step, double *subBoundary,double *minimum_memory,double im_resolution, uint8 iteration,int blunder_selected_level,int Py_combined_level,int pyramid_step,double *minmaxHeight,double ***RPCs,double **ImageAdjust,D2DPOINT *Startpos,D2DPOINT* GridPts, D2DPOINT* Grid_wgs,CSize *Imagesizes_ori,CSize **Imagesizes);

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


double** OpenXMLFile_ortho(char* _filename, double* gsd_r, double* gsd_c);
CSize Envihdr_reader_ortho(char *filename);
CSize Envihdr_reader_DEM_ortho(TransParam param, char *filename, double *minX, double *maxY, double *grid_size);
char* remove_ext_ortho(char* mystr);


//LSF smoothing
CSize GetDEMsize(char *GIMP_path, char* metafilename,TransParam* param, double *grid_size, float* seeddem, double* _minX, double* _maxY);
float* GetDEMValue(char *GIMP_path,CSize seeddem_size);
unsigned char* GetMatchtagValue(char *GIMP_path,CSize seeddem_size);
void DEM_STDKenel_LSF(CSize seeddem_size, bool check_smooth_iter, double MPP_stereo_angle, LSFINFO *Grid_info, double* sigma_average,double* sigma_std, int smooth_iteration,double grid_size,float *seeddem, float *smooth_DEM);
double LocalSurfaceFitting_DEM(double MPP, double sigma_th, int smooth_iter, LSFINFO *Grid_info, float *input, int row_size, int col_size, double grid, long int X, long int Y, long int *numpts, double *fitted_Z);
void LSFSmoothing_DEM(char *savepath, char* outputpath, TransParam param, bool Hemisphere, double MPP, double grid_size, int divide);
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

D2DPOINT GetPhotoCoordinate_single(D3DPOINT A, EO Photo, CAMERA_INFO Camera, RM M);
D3DPOINT GetObjectCoordinate_single(D2DPOINT a, double z,EO Photo, CAMERA_INFO Camera, RM M);
D2DPOINT PhotoToImage_single(D2DPOINT _photo, float _CCDSize, CSize _imgsize);
D2DPOINT ImageToPhoto_single(D2DPOINT _image, float _CCDSize, CSize _imgsize);

bool OpenDMCproject(char* project_path,ProInfo *proinfo, ARGINFO args);
void SetDEMBoundary_photo(EO Photo, CAMERA_INFO m_Camera, RM M, double* _boundary, double* _minmaxheight, double* _Hinterval);
bool SetDEMBoundary_ortho_photo(CSize *Imagesize, double *Boundary, double gridspace, CSize DEM_size, double minX, double maxY, double Ortho_resolution, EO Photo, CAMERA_INFO m_Camera, RM M);
double Correlate(double *L, double *R, int N);
double InterpolatePatch(uint16 *Image, long int position, CSize ImageSize, double dx, double dy);


//Image Coregistration
double** ImageCoregistration(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath, int gcp_opt, D2DPOINT *adjust_std);
void Preprocessing_Coreg(ProInfo *proinfo, char *save_path,uint16 **Oriimage,char **Subsetfile, uint8 py_level, CSize *Subsetsize, CSize **data_size_lr);
int* CoregParam_Image(ProInfo *proinfo,uint8 Pyramid_step, uint8 total_level, double **ImageAdjust, NCCflag _flag,
                      uint8 Template_size, uint16 **Images, CSize **Imagesizes, double **Boundary, double *grid_dx, double *grid_dy,
                      int* grid_space,double** over_Boundary, char* save_path, double* avg_rho, int* iter_count, D2DPOINT *adjust_std);
bool postNCC_ortho(uint8 Pyramid_step, double Ori_diff, double Left_CR,  double Left_CC, double Right_CR, double Right_CC, double **subA,double **TsubA,double **InverseSubA, uint8 Template_size,
                   NCCflag _flag, double bin_angle, CSize leftsize, CSize rightsize, uint16* _leftimage, uint16* _rightimage, double *sum_weight_X, double *sum_weight_Y, double *sum_max_roh);
double *Readtiff_Coreg(char *filename, CSize *Imagesize, int *cols, int *rows, CSize *data_size);
double* LoadPyramidImages_double(char *save_path,char *subsetfile, CSize data_size, uint8 py_level);
double* CreateImagePyramid_double(double* _input, CSize _img_size, int _filter_size, double _sigma);
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
uint16* CreateImagePyramid_avg(uint16* _input, CSize _img_size, int _filter_size);
unsigned char* CreateImagePyramid_BYTE(unsigned char* _input, CSize _img_size, int _filter_size, double _sigma);
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
int Matching_SETSM_SDM(ProInfo proinfo,uint8 pyramid_step, uint8 Template_size, double *Image_res,double *Res, double *Limageparam, double *Rimageparam, CSize Limagesize,CSize Rimagesize, double *Boundary, ImageGSD gsd_image1, ImageGSD gsd_image2,double* LBoundary,double* RBoundary, int *matching_number);
bool subsetImage_SDM(ProInfo proinfo,char *LImageFilename, char *RImageFilename, double *subBoundary, D2DPOINT *Lstartpos, D2DPOINT *Rstartpos, char *LsubsetImage, char *RsubsetImage, CSize* Lsubsetsize, CSize* Rsubsetsize, FILE *fid,bool check_checktiff);
bool GetsubareaImage_SDM(ProInfo proinfo, char *ImageFilename, CSize *Imagesize, double *subBoundary, int *cols, int *rows);
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
void echoprint_Gridinfo_SDM(uint8 prc_level, ProInfo proinfo, double* LBoundary,double* RBoundary, CSize LImagesize, uint16* LeftImage, CSize RImagesize, uint16* RightImage,double* boundary,char *save_path,int row,int col,int level, int iteration, CSize *Size_Grid2D, UGRIDSDM *GridPT3, char *add_str);
void echoprint_adjustXYZ(uint8 prc_level, double* LBoundary,double* RBoundary, CSize LImagesize, uint16* LeftImage, CSize RImagesize, uint16* RightImage,double* boundary,ProInfo proinfo, char *save_path,int row,int col,int level, int iteration, CSize *Size_Grid2D, UGRIDSDM *GridPT3, char *add_str, double gridsize, int d_date);
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
#endif // SETSM_CODE_H
