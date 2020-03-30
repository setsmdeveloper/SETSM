//
//  Coregistration.hpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//

#ifndef Coregistration_hpp
#define Coregistration_hpp

#include "SubFunctions.hpp"

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
F3DPOINT FindNormal(F3DPOINT *normal_ori, float* dem, F3DPOINT Pos,F3DPOINT Mean, F3DPOINT Scale, double* Boundary, Conformalparam X, CSize tinsize, double Gridspace, double *roh_array, float *Z, bool check_tar);
F3DPOINT SurfaceDistance_ori(F3DPOINT tar_normal_ori, float* ref_dem,F3DPOINT tar_normal, F3DPOINT tar_pts, double *tin_boundary,CSize tinsize, double Gridspace, Conformalparam param, F3DPOINT Mean, F3DPOINT Scale, float p_ref_z);
void SlopeAspect(F3DPOINT normal, F3DPOINT scale, float *slope, float *aspect);
F3DPOINT ConformalTransform(F3DPOINT input, Conformalparam param);
F3DPOINT Normalize_coord(F3DPOINT input, F3DPOINT Mean, F3DPOINT Scale);
F3DPOINT Denormalize_coord(F3DPOINT input, F3DPOINT Mean, F3DPOINT Scale);

unsigned char* CreateHillshade(float* _input, CSize _img_size, double grid_size);
F3DPOINT* SettingControls(float* DEM_ref, float* DEM_tar, double grid_size_ref, double grid_size_tar, double *boundary_ref, double *boundary_tar, double* overlapped_br, CSize img_size_ref, CSize img_size_tar, long *tin_point_num);
F3DPOINT* CreateImagePyramid_DEM(float* _input, double grid_size, double *boundary, double* overlapped_br, uint8 pyramid_level, CSize _img_size, int _filter_size, double _sigma, float* result_img, long *tin_point_num, bool check_pts);
void SetHeightRange_slope_aspect(float* ref_img, double* ref_br, CSize ref_size, float* tar_img, double* tar_br, CSize tar_size, long numOfPts, long numOfTri, F3DPOINT *pts, UI3DPOINT *tris, double *boundary, CSize input_size, double input_grid,TINinfo* tininfo, bool check_tar);

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


#endif /* Coregistration_hpp */
