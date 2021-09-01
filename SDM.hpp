//
//  SDM.hpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//

#ifndef SDM_hpp
#define SDM_hpp

#include "SubFunctions.hpp"

//SDM ortho
bool SDM_ortho(char* _filename, ARGINFO args, double** Coreg_param);

void Matching_SETSM_SDM(ProInfo proinfo, TransParam param, uint8 Template_size, double *Rimageparam, const CSize Limagesize, const CSize Rimagesize, const double *Boundary, const ImageGSD gsd_image1, const ImageGSD gsd_image2, long *matching_number, UI2DPOINT &iter_row, UI2DPOINT &iter_col, double buffer_area, double subX, double subY);

bool subsetImage_SDM(ProInfo proinfo, double *subBoundary, D2DPOINT *startpos, CSize* subsetsize, uint16 **Sourceimage);

UGRIDSDM *SetGrid3PT_SDM(const ProInfo *proinfo, LevelInfo &rlevelinfo, const CSize Size_Grid2D, const double Th_roh);

void SetSDMWithSeed(const ProInfo proinfo, LevelInfo &rlevelinfo, UGRIDSDM *Grid, int option);

void SetGridHeightFromSeed(const ProInfo proinfo, LevelInfo &rlevelinfo, UGRIDSDM *Grid, float *seeddem, CSize seeddem_size, double seed_grid, double minX, double maxY, int option);

void SetThs_SDM(int level, int final_level_iteration, double *Th_roh, double *Th_roh_min, double *Th_roh_next, double *Th_roh_start, uint8 pyramid_step);

D2DPOINT *SetGrids_SDM(ProInfo proinfo, const int prc_level, const int level, const int start_py, const int final_level_iteration, CSize *Size_Grid2D, double *py_resolution, double *grid_resolution, const double *subBoundary);

UGRIDSDM *SetGrid3PT_SDM(const ProInfo proinfo, LevelInfo &rlevelinfo, const CSize Size_Grid2D, const double Th_roh);

UGRIDSDM* ResizeGirdPT3_SDM(const CSize preSize, const CSize resize_Size, const double* Boundary, const D2DPOINT *resize_Grid, UGRIDSDM *preGridPT3, const double pre_gridsize);

bool VerticalLineLocus_SDM(ProInfo proinfo, LevelInfo &plevelinfo, NCCresultSDM* nccresult, UGRIDSDM *GridPT3, const ImageGSD gsd_image1, const ImageGSD gsd_image2, double* Coreg_param);

long SelectMPs_SDM(ProInfo proinfo, LevelInfo &rlevelinfo, NCCresultSDM* roh_height, UGRIDSDM *GridPT3, vector<D3DPOINT> &Matched_pts_col, vector<D3DPOINT> &Matched_pts_row);

void echoprint_Gridinfo_SDM(ProInfo proinfo, LevelInfo &rlevelinfo, int row, int col, int level, int iteration, UGRIDSDM *GridPT3);

void echoprint_adjustXYZ(ProInfo proinfo, LevelInfo &rlevelinfo, int row, int col,int level, int iteration, UGRIDSDM *GridPT3, int d_date);

bool Update_ortho_NCC(ProInfo proinfo, LevelInfo &rlevelinfo, UGRIDSDM *GridPT3, const ImageGSD gsd_image1, const ImageGSD gsd_image2, double* Coreg_param);

UGRIDSDM* CopyGridPT3(UGRIDSDM *GridPT3, LevelInfo &rlevelinfo);

void SetShiftFromTIN(const long numOfPts, const long numOfTri, UGRIDSDM *GridPT3, LevelInfo &rlevelinfo, D3DPOINT *pts, UI3DPOINT *tris, const int b_dir);

void shift_filtering(ProInfo proinfo, UGRIDSDM *GridPT3, LevelInfo &rlevelinfo);

//unused functions
void echo_print_nccresults_SDM(char *save_path,int row,int col,int level, int iteration, NCCresultSDM *nccresult, CSize *Size_Grid2D, char *add_str);

bool average_filter_colrowshift(CSize Size_Grid2D, UGRIDSDM *GridPT3,uint8 Pyramid_step);

double MergeTiles_SDM(ProInfo info,int iter_row_end,int t_col_end, double buffer,int final_iteration,TransParam _param, uint8 pyramid_step);

#endif /* SDM_hpp */
