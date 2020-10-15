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


#include <omp.h>
#include <time.h>
#include <syslog.h>
#include <unistd.h>

#include "git_description.h"

#include "SubFunctions.hpp"
#include "LSF.hpp"
#include "Orthogeneration.hpp"
#include "Coregistration.hpp"
#include "SDM.hpp"


void DownSample(ARGINFO &args);

int SETSMmainfunction(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath,double **ImageParam);

bool SetupParam(ProInfo *info, bool *pre_DEMtif);

void SetTiles(const ProInfo *info, const double *Boundary, const int tile_size, uint8 *pyramid_step, uint16 *buffer_area, uint8 *iter_row_start, uint8 *iter_row_end, uint8 *t_col_start, uint8 *t_col_end, double *subX, double *subY);

void SetTiles_RA(const ProInfo *info, const double *Boundary, const int tile_size, uint8 *pyramid_step, uint8 *RA_row_start, uint8 *RA_row_end, uint8 * RA_row_iter, uint8 *t_col_start, uint8 *t_col_end, uint8 *RA_col_iter, double *subX, double *subY);

void SetsubsetBR(const ProInfo *info, const double *Boundary, const int tile_size, double *subX, double *subY, int &division_X, int &division_Y);

void SetThs(const ProInfo *proinfo,const int level, const int final_level_iteration, double *Th_roh, double *Th_roh_min,double *Th_roh_next, double *Th_roh_start);

D2DPOINT *SetGrids(const ProInfo *info, const int level, const int final_level_iteration, const double resolution, CSize *Size_Grid2D, const double DEM_resolution, double *py_resolution, double *grid_resolution, const double *subBoundary);

UGRID *SetGrid3PT(const ProInfo *proinfo, LevelInfo &rlevelinfo, const double Th_roh, double *minmaxHeight);

void SetSubBoundary(const double *Boundary, const double subX, const double subY, const double buffer_area, const int col, const int row, double *subBoundary);

void SetHeightWithSeedDEM(const ProInfo *proinfo, LevelInfo &rlevelinfo, UGRID *Grid, double *minmaxHeight);

void SetGridHeightFromSeed(LevelInfo &rlevelinfo, UGRID *Grid, float *seeddem, CSize seeddem_size, double seed_grid, double minX, double maxY, double seedDEM_sigma, double *minmaxHeight);

void SetDEMBoundary(double** _rpcs, double* _res,TransParam _param, double* _boundary, double* _minmaxheight, double* _Hinterval);

void SetDEMBoundary_photo(EO Photo, CAMERA_INFO m_Camera, RM M, double* _boundary, double* _minmaxheight, double* _Hinterval);

uint16 *SetsubsetImage(ProInfo *proinfo, LevelInfo &rlevelinfo, const int index_image, const TransParam transparam, const uint8 NumofIAparam, const double * const * const *RPCs, const double * const *ImageParams, const double *subBoundary, const double *minmaxHeight, D2DPOINT *Startpos, CSize *Subsetsize);

void CalMPP_pair(double CA,double mean_product_res, double im_resolution, double *MPP_stereo_angle);

void CalMPP(ProInfo *proinfo, LevelInfo &rlevelinfo, const double* minmaxHeight, const double CA,const double mean_product_res, double *MPP_simgle_image, double *MPP_stereo_angle);

void CalMPP_8(ProInfo *proinfo, LevelInfo &rlevelinfo, const double* minmaxHeight, const double CA,const double mean_product_res, double *MPP_simgle_image, double *MPP_stereo_angle);

void InitializeVoxel(const ProInfo *proinfo, VOXEL **grid_voxel,LevelInfo &plevelinfo, UGRID *GridPT3, NCCresult* nccresult,const int iteration, const double *minmaxHeight);

double GetHeightStep(int Pyramid_step, double im_resolution);

double SetNCC_alpha(const int Pyramid_step, const int iteration, bool IsRA);

int VerticalLineLocus(VOXEL **grid_voxel,const ProInfo *proinfo, NCCresult* nccresult, LevelInfo &plevelinfo, const UGRID *GridPT3, const uint8 iteration,const double *minmaxHeight);

void SetOrthoImageCoord(const ProInfo *proinfo, LevelInfo &plevelinfo, const UGRID *GridPT3, const bool check_combined_WNCC, enum PyImageSelect check_pyimage, const double im_resolution, const double im_resolution_next, long int &sub_imagesize_w, long int &sub_imagesize_h, long int &sub_imagesize_w_next, long int &sub_imagesize_h_next, D2DPOINT **am_im_cd, D2DPOINT **am_im_cd_next);

void FindPeakNcc(const int Pyramid_step, const int iteration, const long int grid_index, const double temp_rho, const float iter_height, bool &check_rho, double &pre_rho, float &pre_height, int &direction, double &max_WNCC, NCCresult *nccresult);

void SGM_start_pos(NCCresult *nccresult, VOXEL** grid_voxel, UGRID *GridPT3, long pt_index, float* LHcost_pre, float **SumCost, double height_step_interval);

void SGM_con_pos(int pts_col, int pts_row, CSize Size_Grid2D, int direction_iter, double step_height, int P_HS_step, int *u_col, int *v_row, NCCresult *nccresult, VOXEL** grid_voxel,UGRID *GridPT3, long pt_index, double P1, double P2, float* LHcost_pre, float* LHcost_curr, float **SumCost);

void AWNCC(ProInfo *proinfo, VOXEL **grid_voxel,CSize Size_Grid2D, UGRID *GridPT3, NCCresult *nccresult, double step_height, uint8 Pyramid_step, uint8 iteration,int MaxNumberofHeightVoxel);

void VerticalLineLocus_seeddem(const ProInfo *proinfo,LevelInfo &rlevelinfo, UGRID *GridPT3, const double* minmaxHeight);

bool VerticalLineLocus_blunder(const ProInfo *proinfo,LevelInfo &rlevelinfo, float* nccresult, UGRID *GridPT3, uint8 iteration, bool bblunder);

int VerticalLineLocus_Ortho(ProInfo *proinfo, LevelInfo &rlevelinfo, double MPP, double *F_Height, D3DPOINT ref1_pt, D3DPOINT ref2_pt, D3DPOINT target_pt, UGRID *GridPT3, int target_index, double *F_sncc);

long SelectMPs(const ProInfo *proinfo,LevelInfo &rlevelinfo, const NCCresult* roh_height, UGRID *GridPT3, const double Th_roh, const double Th_roh_min, const double Th_roh_start, const double Th_roh_next, const int iteration, const double MPP, const int final_level_iteration,const double MPP_stereo_angle, vector<D3DPOINT> *linkedlist);

UI3DPOINT* TINgeneration(bool last_flag, char *savepath, uint8 level, CSize Size_Grid2D, double img_resolution, double grid_resolution,
						 double min_max[],
						 double *subBoundary, int total_point_count, D3DPOINT *ptslists, int *iter_row, int *iter_col,
						 int *re_total_tri_counts);

void DecisionMPs(const ProInfo *proinfo, LevelInfo &rlevelinfo, const bool flag_blunder,const long int count_MPs_input,UGRID *GridPT3, const uint8 iteration, const double Hinterval, int *count_Results, double *minz_mp, double *maxz_mp, const double *minmaxHeight, D3DPOINT *ptslists);

void DecisionMPs_setheight(const ProInfo *proinfo, LevelInfo &rlevelinfo, const long int count_MPs_input,UGRID *GridPT3, const uint8 iteration, const double Hinterval, const double *minmaxHeight, D3DPOINT *ptslists, UI3DPOINT *trilists,int numoftri);

int SetttingFlagOfGrid(LevelInfo &rlevelinfo, UGRID *GridPT3, vector<D3DPOINT> MatchedPts_list_anchor,vector<D3DPOINT> MatchedPts_list_blunder,vector<D3DPOINT> *MatchedPts_list_mps);

int AdjustParam(ProInfo *proinfo, LevelInfo &rlevelinfo, int NumofPts, double **ImageAdjust, uint8 total_pyramid, D3DPOINT* ptslists);

bool postNCC(LevelInfo &rlevelinfo, const double Ori_diff, const D2DPOINT left_pt, const D2DPOINT right_pt, double subA[][6], double TsubA[][9], double InverseSubA[][6], uint8 Half_template_size, const int reference_ID, const int target_ID, double *sum_weight_X, double *sum_weight_Y, double *sum_max_roh, double **left_patch_vecs, double **right_patch_vecs);

bool blunder_detection_TIN(const ProInfo *proinfo, LevelInfo &rlevelinfo, const int iteration, float* ortho_ncc, bool flag_blunder, uint16 count_bl, D3DPOINT *pts, bool *detectedBlunders, long int num_points, UI3DPOINT *tris, long int num_triangles, UGRID *Gridpts, long *blunder_count,double *minz_mp, double *maxz_mp);

int Ortho_blunder(ProInfo *proinfo, LevelInfo &rlevelinfo, double MPP, D3DPOINT *pts, int numOfPts, UI3DPOINT *tris,int numOfTri, UGRID *GridPT3);

bool SetHeightRange_blunder(LevelInfo &rlevelinfo, const D3DPOINT *pts, const int numPts, UI3DPOINT *tris,const long num_triangles, UGRID *GridPT3);

UGRID* SetHeightRange(ProInfo *proinfo, LevelInfo &rlevelinfo, NCCresult *nccresult, const int numOfPts, const int num_triangles, UGRID *GridPT3, const int iteration, double *minH_grid, double *maxH_grid, D3DPOINT *pts, const UI3DPOINT *tris, const double MPP, const bool level_check_matching_rate);

UGRID* ResizeGirdPT3(ProInfo *proinfo, CSize preSize, CSize resize_Size, double* Boundary, D2DPOINT *resize_Grid, UGRID *preGridPT3, double pre_gridsize, double* minmaxheight);
UGRID* ResizeGirdPT3_RA(const ProInfo *proinfo,const CSize preSize,const CSize resize_Size,const double* preBoundary,const double* Boundary,const D2DPOINT *resize_Grid, UGRID *preGridPT3,const double pre_gridsize,const double* minmaxheight);

void echoprint_Gridinfo(ProInfo *proinfo, NCCresult* roh_height, int row,int col,int level, int iteration, double update_flag, CSize *Size_Grid2D, UGRID *GridPT3, char *add_str);
void echo_print_nccresults(char *save_path,int row,int col,int level, int iteration, NCCresult *nccresult, CSize *Size_Grid2D, char *add_str);

int Matching_SETSM(ProInfo *proinfo,const uint8 pyramid_step, const uint8 Template_size, const uint16 buffer_area, const uint8 iter_row_start, const uint8 iter_row_end, const uint8 t_col_start, const uint8 t_col_end, const double subX,const double subY,const double bin_angle,const double Hinterval,const double *Image_res, double **Imageparams, const double *const*const*RPCs, const uint8 NumOfIAparam, const CSize *Imagesizes,const TransParam param, double *minmaxHeight,const double *Boundary, const double CA,const double mean_product_res, double *stereo_angle_accuracy);

bool check_kernel_size(ProInfo *proinfo, const CSize *Subsetsize,const int Template_size, const int pyramid_step);

bool check_image_boundary(const ProInfo *proinfo,LevelInfo &plevelinfo, const D2DPOINT pos_xy_m,const D2DPOINT pos_xy,const double minH,const double maxH,const int H_template_size);

double CalMemorySize_Post(CSize DEM_size,CSize Final_DEMsize);

double CalMemorySize_Post_MT(CSize DEM_size, CSize Final_DEMsize);

double CalMemorySize_Post_LSF(CSize DEM_size, CSize Final_DEMsize);

void MergeTiles(const ProInfo *info, const int iter_row_start, const int t_col_start, const int iter_row_end, const int t_col_end, int buffer, const int final_iteration, float *DEM, const CSize Final_DEMsize, double *FinalDEM_boundary);

void MergeTiles_Ortho(const ProInfo *info, const int iter_row_start, const int t_col_start, const int iter_row_end,const int t_col_end, int buffer,const int final_iteration, signed char *DEM_ortho, const CSize Final_DEMsize, const double *FinalDEM_boundary);

double FindNebPts_F_M_IDW(const float *input, const unsigned char *matching_flag, const long row_size, const long col_size, const double grid, const double minX, const double maxY, const double X, const double Y, const int row_interval);

CSize DEM_final_Size(const char *save_path, const int row_start, const int col_start,const int row_end, const int col_end, const double grid_resolution, double *boundary);

void NNA_M(const ProInfo *proinfo, const TransParam _param, const int row_start, const int col_start, const int row_end, const int col_end, int buffer_clip, const int final_iteration, const int divide, const CSize Final_DEMsize, float* DEM_values, float* value, unsigned char* value_pt, const double *FinalDEM_boundary);

void NNA_M_MT(const ProInfo *proinfo, const TransParam _param, const int row_start, const int col_start,const int row_end, const int col_end, int buffer_clip, const int final_iteration, const int divide, signed char* Ortho_values, float* value, unsigned char* value_pt, const CSize Final_DEMsize, const double *FinalDEM_boundary);

double CalMemorySize(const ProInfo *info,LevelInfo &plevelinfo,const UGRID *GridPT3, double *minimum_memory, const uint8 iteration, const double *minmaxHeight);

class TileIndexer {
    public:
    virtual int next() = 0;
    virtual ~TileIndexer() {};
};

struct SINGLE_FILE {
    FILE *fp;
};

SINGLE_FILE *single_fopen(const char *path, const char *mode);
int fclose(SINGLE_FILE *stream);
int fprintf(SINGLE_FILE *stream, const char *format, ...);
int single_printf(const char *fmt, ...);


#endif // SETSM_CODE_H
