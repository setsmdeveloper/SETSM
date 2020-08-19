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
//End Image Coregistration

void DEM_ImageCoregistration_hillshade(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath, int gcp_opt);
void DEM_ImageCoregistration_GeomatricConstraint(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath, int gcp_opt);
D3DPOINT FindNormal(D3DPOINT *normal_ori, const float* dem, const D3DPOINT Pos, const D3DPOINT Mean, const D3DPOINT Scale, const double* Boundary, const Conformalparam X, const CSize tinsize, const double Gridspace, double *roh_array, double *Z, const bool check_tar);
D3DPOINT SurfaceDistance_ori(const D3DPOINT tar_normal_ori, const float* ref_dem, const D3DPOINT tar_normal, const D3DPOINT tar_pts, const double *tin_boundary,const CSize tinsize, const double Gridspace, const Conformalparam param, const D3DPOINT Mean, const D3DPOINT Scale, const double p_ref_z);
void SlopeAspect(const D3DPOINT normal, const D3DPOINT scale, double *slope, double *aspect);
D3DPOINT ConformalTransform(const D3DPOINT input, const Conformalparam param);
D3DPOINT Normalize_coord(const D3DPOINT input, const D3DPOINT Mean, const D3DPOINT Scale);
D3DPOINT Denormalize_coord(const D3DPOINT input, const D3DPOINT Mean, const D3DPOINT Scale);

unsigned char* CreateHillshade(const float* _input, const CSize _img_size, const double grid_size);
void SettingControls(const ProInfo *proinfo, const float* DEM_ref, const float* DEM_tar, const double grid_size_ref, const double grid_size_tar, const double *boundary_ref, const double *boundary_tar, const double* overlapped_br, const CSize img_size_ref, const CSize img_size_tar, long &tin_point_num, vector<D3DPOINT> &select_pts);

double* CoeffMatrix_25D(const D3DPOINT coord_center, const D3DPOINT coord_scale, const long selected_pts, const vector<D3DPOINT> &transformed_coord, const vector<double> &dH, const Conformalparam param, const vector<D3DPOINT> &tar_normal, const vector<double> &weight, double* sigmaX, double *sigma0);


#endif /* Coregistration_hpp */
