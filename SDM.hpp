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

void Matching_SETSM_SDM(const ProInfo proinfo, uint8 Template_size, double *Rimageparam, const CSize Limagesize, const CSize Rimagesize, const double *Boundary, const ImageGSD gsd_image1, const ImageGSD gsd_image2, int *matching_number);

void SetHeightWithSeedDEM_SDM(ProInfo proinfo,UGRIDSDM *Grid, double *Boundary, CSize Grid_size, double Grid_set);

bool subsetImage_SDM(ProInfo proinfo, double *subBoundary, D2DPOINT *startpos, CSize* subsetsize, uint16 **Sourceimage);

bool GetsubareaImage_SDM(ProInfo proinfo, char *ImageFilename, CSize *Imagesize, double *subBoundary, long int *cols,long int *rows);
D2DPOINT* GetObjectToImage(uint16 _numofpts, D2DPOINT *_GP, double *boundary, double imageres);
void Preprocessing_SDM(ProInfo proinfo, char *save_path,char *Lsubsetfile, char *Rsubsetfile, uint8 py_level, CSize *Lsubsetsize, CSize *Rsubsetsize, CSize *data_size_l, CSize *data_size_r, FILE *fid);
void SetThs_SDM(int level, int final_level_iteration, double *Th_roh, double *Th_roh_min, double *Th_roh_next, double *Th_roh_start, uint8 pyramid_step);

D2DPOINT *SetGrids_SDM(ProInfo proinfo, const int prc_level, const int level, const int start_py, const int final_level_iteration, CSize *Size_Grid2D, double *py_resolution, double *grid_resolution, const double *subBoundary);

UGRIDSDM *SetGrid3PT_SDM(const CSize Size_Grid2D, const double Th_roh);

UGRIDSDM* ResizeGirdPT3_SDM(const CSize preSize, const CSize resize_Size, const double* Boundary, const D2DPOINT *resize_Grid, UGRIDSDM *preGridPT3, const double pre_gridsize);

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
double MergeTiles_SDM(ProInfo info,int iter_row_end,int t_col_end, int buffer,int final_iteration,TransParam _param, uint8 pyramid_step);





#endif /* SDM_hpp */
