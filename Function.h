#include "Typedefine.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "tiff.h"
#include "tiffio.h"
//#include "vdefs.h"

#define  MAXRAND     0x7fffffff
#define  BIGNUM      1e37
#define  WEEBIT      0.000000000001
#define  MAXDIM      10
#define  MAXSTR      48
#define  SQ(x)       (x) * (x)

FILE *fid_bisector, *fid_ep, *fid_vertex, *fid_site, *fid_triple;
int sorted, triangulate, plot, debug, count_tri;
int nsites, siteidx ;
float xmin, xmax, ymin, ymax ;
Site * sites ;
Freelist sfl ;

struct simplex
{  int            *vert;
   double         *cent;
   struct simplex *nextsimp;
};

struct facelist
{  int             *vert;
   struct facelist *nextface;
};

struct nnweight       
{  double          weight;
   int             vert;
   struct nnweight *nextcoord;
};

typedef struct nnXY
{
    double X;
    double Y;
    double Z;
} NNXY;

void SETSMmainfunction(char* _filename, ARGINFO args, char *_LeftImagefilename, char *_save_filepath);
bool OpenProject_all(char* _filename, ProInfo *info);
bool OpenProject(char* _filename, ProInfo *info, ARGINFO args);
int Maketmpfolders(ProInfo info);
bool SetupParam_all(ProInfo info,uint8 *NumOfIAparam, uint8 *pre_DEM_level, uint8 *DEM_level,  bool *pre_DEMtif, bool *check_tile_array, bool *tile_array, float *Boundary);
bool SetupParam(ProInfo info,uint8 *NumOfIAparam, uint8 *pre_DEM_level, uint8 *DEM_level,  bool *pre_DEMtif, bool *check_tile_array);
void SetTransParam(float minLat, bool *Hemisphere, TransParam *param);
void SetTiles(ProInfo info, bool IsSP, bool IsRR, float *Boundary, float *Res, int tile_size, bool pre_DEMtif, uint8 *pyramid_step, uint16 *buffer_area, 
			  uint8 *iter_row_start, uint8 *iter_row_end, uint8 *t_col_start, uint8 *t_col_end, float *subX, float *subY);
void SetTiles_RA(ProInfo info, bool IsSP, bool IsRR, float *Boundary, float *Res, int tile_size, bool pre_DEMtif, uint8 *pyramid_step, uint16 *buffer_area, 
			  uint8 *RA_row_start, uint8 *RA_row_end, uint8 * RA_row_iter, uint8 *t_col_start, uint8 *t_col_end, uint8 *RA_col_iter, float *subX, float *subY);
void SetPySizes(CSize *data_size_l, CSize *data_size_r, CSize Lsubsetsize, CSize Rsubsetsize, int level);
void SetThs_ratio(int level, float *Th_roh, float *Th_roh_min, float *Th_roh_next, float *Th_roh_start, int pre_DEMtif, int IsRA);
void SetThs(int level, float *Th_roh, float *Th_roh_min, float *Th_roh_next, float *Th_roh_start, int pre_DEMtif, int IsRA, float seedDEMsigma);
F2DPOINT *SetGrids(bool *dem_update_flag, bool flag_start, int level, float resolution, CSize *Size_Grid2D, bool pre_DEMtif, char *priori_DEM_tif, float DEM_resolution, float *minmaxHeight,
				   float *py_resolution, float *grid_resolution, float *subBoundary);
UGRID *SetGrid3PT(bool dem_update_flag, bool flag_start, CSize Size_Grid2D, float Th_roh, int level, float *minmaxHeight,float *subBoundary,float py_resolution,
				  char* priori_DEM_tif,bool pre_DEMtif, float seedDEMsigma, int IsRA,char* metafilename);
int	 CalTotalIteration(uint8 DEM_level,int level);

char* remove_ext(char* mystr);

char* GetFileName(char file_path[]);
char* GetFileDir(char file_path[],int *size);

bool GetImageSize(char *filename, CSize *Imagesize);
bool GetsubareaImage(TransParam transparam, uint8 NumofIAparam, double **RPCs, float *ImageParam, char *ImageFilename, CSize *Imagesize,
				 float *subBoundary, float *minmaxHeight, int *cols, int *rows);
uint16 *Readtiff(char *filename, CSize *Imagesize, int *cols, int *rows, CSize *data_size,bool check_checktiff);
bool Writetiff(char *filename, float* input, CSize Imagesize);

float *Readtiff_DEM(char *filename, CSize *Imagesize, int *cols, int *rows, CSize *data_size);
void SetSubBoundary(float *Boundary, float subX, float subY, float buffer_area, int col, int row, float *subBoundary);
F2DPOINT *SetDEMGrid(float *Boundary, float Grid_x, float Grid_y, CSize *Size_2D);
void SetHeightWithSeedDEM(UGRID *Grid, float *Boundary,CSize Grid_size, float Grid_set,  char *GIMP_path, float *minmaxHeight, float seedDEMsigma, int IsRA,char* metafilename);
double** OpenXMLFile(char* _filename,double* gsd_r, double* gsd_c);
void SetDEMBoundary(double** _rpcs, float* _res,TransParam _param, bool _hemisphere, float* _boundary, float* _minmaxheight, CSize* _imagesize, float* _Hinterval);
bool subsetImage(TransParam transparam, uint8 NumofIAparam, double **LRPCs, float *LImageParam, char *LImageFilename, double **RRPCs, float *RImageParam, char *RImageFilename, 
				 float *subBoundary, float *minmaxHeight, F2DPOINT *Lstartpos, F2DPOINT *Rstartpos, char *LsubsetImage, char *RsubsetImage, CSize* Lsubsetsize, CSize* Rsubsetsize, FILE *fid,bool check_checktiff);
F2DPOINT* wgs2ps(bool _bopenmp, TransParam _param, int _numofpts, F2DPOINT *_wgs);
F2DPOINT wgs2ps_single(TransParam _param, F2DPOINT _wgs);
F3DPOINT* wgs2ps_3D(bool _bopenmp, TransParam _param, int _numofpts, F3DPOINT *_wgs);
F2DPOINT* ps2wgs(bool _bopenmp, TransParam _param, int _numofpts, F2DPOINT *_ps);
F2DPOINT ps2wgs_single(TransParam _param, F2DPOINT _ps);
F3DPOINT* ps2wgs_3D(bool _bopenmp, TransParam _param, int _numofpts, F3DPOINT *_ps);
F2DPOINT* GetObjectToImageRPC(double **_rpc, uint8 _numofparam, float *_imageparam, uint16 _numofpts, F3DPOINT *_GP);
F2DPOINT GetObjectToImageRPC_single(double **_rpc, uint8 _numofparam, float *_imageparam, F3DPOINT _GP);
void Preprocessing(char *save_path,char *Lsubsetfile, char *Rsubsetfile, uint8 py_level, CSize *Lsubsetsize, CSize *Rsubsetsize, CSize *data_size_l, CSize *data_size_r, FILE *fid);
uint16* LoadPyramidImages(char *save_path,char *subsetfile, CSize data_size, uint8 py_level);
uint16* LoadPyramidMagImages(char *save_path,char *subsetfile, CSize data_size, uint8 py_level);
uint8* LoadPyramidOriImages(char *save_path,char *subsetfile, CSize data_size, uint8 py_level);
uint16* CreateImagePyramid(uint16* _input, CSize _img_size, int _filter_size, float _sigma);
void MakeSobelMagnitudeImage(CSize _img_size, uint16* _src_image, uint16* _dist_mag_image, /*int16* _gx, int16* _gy,*/ int16* _dir);
void Orientation(CSize imagesize, uint16* Gmag, int16* Gdir, uint8 Template_size, uint8* plhs);
NCCresult* VerticalLineLocus(uint16 *MagImages_L,uint16 *MagImages_R,float DEM_resolution, float im_resolution, double** LRPCs, double** RRPCs, CSize LImagesize_ori, CSize LImagesize, uint16* LeftImage, CSize RImagesize_ori, CSize RImagesize, uint16* RightImage, uint8 Template_size,
                             CSize Size_Grid2D, TransParam param, F2DPOINT* GridPts, F2DPOINT* Grid_wgs, UGRID *GridPT3, NCCflag flag,
                             uint8 NumofIAparam, float* ImageAdjust, float* minmaxHeight, uint8 Pyramid_step, F2DPOINT Lstartpos, F2DPOINT Rstartpos, uint8 iteration, uint8* left_ori, uint8* right_ori,
                             float bin_angle, uint8 NumOfCompute, uint8 peak_level, FILE* fid, bool IsPar, bool Hemisphere, char* save_filepath, uint8 tile_row, uint8 tile_col, float* Boundary,
                             bool pre_DEMtif, char* tmpdir,float *meters_per_pixel, bool IsRA);

float VerticalLineLocus_seeddem(uint16 *MagImages_L,uint16 *MagImages_R,float DEM_resolution, float im_resolution, double** LRPCs, double** RRPCs,
							   CSize LImagesize_ori, CSize LImagesize, uint16* LeftImage, CSize RImagesize_ori, CSize RImagesize, uint16* RightImage, uint8 Template_size, 
							   CSize Size_Grid2D, TransParam param, F2DPOINT* GridPts, F2DPOINT *Grid_wgs, UGRID *GridPT3,
							   uint8 NumofIAparam, float* ImageAdjust, uint8 Pyramid_step, F2DPOINT Lstartpos, F2DPOINT Rstartpos, 
							   char* save_filepath, uint8 tile_row, uint8 tile_col, uint8 iteration,uint8 bl_count,float* Boundary, float* minmaxHeight, float seedDEMsigma);

bool VerticalLineLocus_blunder(float* nccresult, uint16 *MagImages_L,uint16 *MagImages_R,float DEM_resolution, float im_resolution, double** LRPCs, double** RRPCs,
                               CSize LImagesize_ori, CSize LImagesize, uint16* LeftImage, CSize RImagesize_ori, CSize RImagesize, uint16* RightImage, uint8 Template_size,
                               CSize Size_Grid2D, TransParam param, F2DPOINT* GridPts, F2DPOINT *Grid_wgs, UGRID *GridPT3,
                               uint8 NumofIAparam, float* ImageAdjust, uint8 Pyramid_step, F2DPOINT Lstartpos, F2DPOINT Rstartpos,
                               char* save_filepath, uint8 tile_row, uint8 tile_col, uint8 iteration,uint8 bl_count,float* Boundary);
bool VerticalLineLocus_Ortho(float *F_height,F3DPOINT ref1_pt, F3DPOINT ref2_pt, F3DPOINT target_pt,
							  uint16 *MagImages_L,uint16 *MagImages_R, uint16* LeftImage, uint16* RightImage,
							  float DEM_resolution, float im_resolution, double** LRPCs, double** RRPCs, 
							  CSize LImagesize,  CSize RImagesize, CSize Size_Grid2D, TransParam param, uint8 NumofIAparam, 
							  float* ImageAdjust, float* minmaxHeight, uint8 Pyramid_step, float meters_per_pixel,
							  F2DPOINT Lstartpos, F2DPOINT Rstartpos, uint8 iteration,  UGRID *GridPT3, int target_index,int ref1_index, int ref2_index,
							  float* boundary);
F2DPOINT* OriginalToPyramid(uint16 numofpts,F2DPOINT* InCoord, F2DPOINT Startpos, uint8 Pyramid_step);
F2DPOINT OriginalToPyramid_single(F2DPOINT InCoord, F2DPOINT Startpos, uint8 Pyramid_step);
F2DPOINT* PyramidToOriginal(uint16 numofpts,F2DPOINT* InCoord, F2DPOINT Startpos, uint8 Pyramid_step);
int SelectMPs(NCCresult* roh_height, CSize Size_Grid2D, F2DPOINT *GridPts_XY, UGRID *GridPT3,
				 float Th_roh, float Th_roh_min, float Th_roh_start, float Th_roh_next, uint8 Pyramid_step, uint8 total_pyramid,
				 uint8 iteration, uint8 peak_level, char *filename_mps, int pre_DEMtif, int IsRA, float MPP);

bool TINgeneration(bool last_flag, char *savepath, uint8 level, CSize Size_Grid2D, float img_resolution, float grid_resolution, 
				   F3DPOINT *scaled_ptslists,
				   float *subBoundary, int total_point_count, F3DPOINT *ptslists, int *iter_row, int *iter_col, 
				   int *re_total_tri_counts);

UI3DPOINT *LoadTrilistFromFile(char *savepath, int interval_row, int interval_col, int total_tri_counts);

int DecisionMPs(bool flag_blunder,int count_MPs, float* Boundary, UGRID *GridPT3, uint8 Pyramid_step, float grid_resolution, 
				 uint8 iteration, CSize Size_Grid2D, char *filename_mps_pre, char *filename_mps, char *filename_tri, float Hinterval, 
				 bool *p_flag, float *pre_3sigma, float *pre_mean, int *count_Results, float *minz_mp, float *maxz_mp, float *minmaxHeight,
				uint16 *MagImages_L,uint16 *MagImages_R,float DEM_resolution, float im_resolution, double** LRPCs, double** RRPCs,
				CSize LImagesize_ori, CSize LImagesize, uint16* LeftImage, CSize RImagesize_ori, CSize RImagesize, uint16* RightImage, uint8 Template_size, 
				TransParam param, F2DPOINT* Grid_wgs,F2DPOINT* GridPts,
				uint8 NumofIAparam, float* ImageAdjust, F2DPOINT Lstartpos, F2DPOINT Rstartpos, 
				char* save_filepath, uint8 tile_row, uint8 tile_col, int pre_DEMtif, int IsRA);
int DecisionMPs_setheight(bool flag_blunder, int count_MPs_input, float* Boundary, UGRID *GridPT3, uint8 Pyramid_step, float grid_resolution, 
						  uint8 iteration, CSize Size_Grid2D, char *filename_mps_pre, char *filename_tri, float Hinterval, 
						  bool *p_flag, float *pre_3sigma, float *pre_mean, int *count_Results, float *minz_mp, float *maxz_mp, float *minmaxHeight,
						  uint16 *MagImages_L,uint16 *MagImages_R,float DEM_resolution, float im_resolution, double** LRPCs, double** RRPCs,
						  CSize LImagesize_ori, CSize LImagesize, uint16* LeftImage, CSize RImagesize_ori, CSize RImagesize, uint16* RightImage, uint8 Template_size, 
						  TransParam param, F2DPOINT* Grid_wgs,F2DPOINT* GridPts,
						  uint8 NumofIAparam, float* ImageAdjust, F2DPOINT Lstartpos, F2DPOINT Rstartpos, 
						  char* save_filepath, uint8 tile_row, uint8 tile_col,F3DPOINT *ptslists, UI3DPOINT *trilists,int numoftri);
int SetttingFlagOfGrid(float *subBoundary,UGRID *GridPT3, uint8 Pyramid_step,float grid_resolution,uint8 iteration,
						CSize Size_Grid2D,char *filename_mps_anchor,char *filename_mps_aft,int count_results_anchor,int count_results_blunder, char *filename_mps);

int AdjustParam(uint8 Pyramid_step, int NumofPts, char * file_pts, F2DPOINT Lstartpos, F2DPOINT Rstartpos, double **LRPCs, double **RRPCs, float *ImageAdjust, NCCflag _flag,
					 uint8 Template_size, uint16 *LeftImage, CSize LImagesize, uint16 *RightImage, CSize RImagesize, uint8 *left_ori, uint8 *right_ori, TransParam param,
					 float bin_angle, uint8 total_pyramid, bool Hemisphere, char* save_filepath, char* tmpdir);
bool postNCC(uint8 Pyramid_step, double Ori_diff, double Left_CR,  double Left_CC, double Right_CR, double Right_CC, float **subA,float **TsubA,float **InverseSubA, uint8 Template_size, 
			NCCflag _flag, float bin_angle, CSize leftsize, CSize rightsize, uint16* _leftimage, uint16* _rightimage, double *sum_weight_X, double *sum_weight_Y, double *sum_max_roh);

int scomp(const void * vs1, const void * vs2);
Site *nextone(void);
void readsites(F3DPOINT *ptslists,int numofpts);
Site *readone(void);
void TINCreate(F3DPOINT *ptslists,char *filename_tri,int numofpts);
bool blunder_detection_TIN(int pre_DEMtif, float* ortho_ncc,bool flag_blunder,uint16 count_bl,float* blunder_dh,char *file_pts, 
						   F3DPOINT *ptslists, int numOfPts, UI3DPOINT *trilists,int numOfTri, UGRID *Gridpts, BL BL_param, 
						   uint32 *blunder_count,float *minz_mp, float *maxz_mp, float *minmaxHeight, int IsRA);

int Ortho_blunder(F3DPOINT *pts, int numOfPts, UI3DPOINT *tris,int numOfTri, bool update_flag,float *minH_grid, float *maxH_grid, BL BL_param,
				  uint16 *MagImages_L,uint16 *MagImages_R, uint16* LeftImage, uint16* RightImage,
				  float DEM_resolution, float im_resolution, double** LRPCs, double** RRPCs, 
				  CSize LImagesize,  CSize RImagesize, CSize Size_Grid2D, TransParam param, uint8 NumofIAparam, 
				  float* ImageAdjust, float* minmaxHeight, uint8 Pyramid_step, float meters_per_pixel,
				  F2DPOINT Lstartpos, F2DPOINT Rstartpos, uint8 iteration,  UGRID *GridPT3, char *filename_mps, char *save_path);

bool SetHeightRange_blunder(float* minmaxHeight,F3DPOINT *pts, int numOfPts, UI3DPOINT *tris,int numOfTri, UGRID *GridPT3, BL BL_param, float *mt_minmaxheight);
UGRID* SetHeightRange(bool pre_DEMtif, float* minmaxHeight,int numOfPts, int numOfTri, UGRID *GridPT3, bool update_flag,
					  float *minH_grid, float *maxH_grid, BL BL_param,F3DPOINT *pts, UI3DPOINT *tris,int IsRA, float MPP);

void echoprint_Gridinfo(char *save_path,int row,int col,int level, int iteration, float update_flag, CSize *Size_Grid2D, UGRID *GridPT3, char *add_str);
void echo_print_nccresults(char *save_path,int row,int col,int level, int iteration, NCCresult *nccresult, CSize *Size_Grid2D, char *add_str);

bool Matching_SETSM(ProInfo proinfo,uint8 pyramid_step, uint8 Template_size, uint16 buffer_area,uint8 iter_row_start, uint8 iter_row_end,uint8 t_col_start,uint8 t_col_end,
			  float subX,float subY,float bin_angle,float Hinterval,float *Image_res,float *Res, float *Limageparam, float *Rimageparam,
			  double **LRPCs, double **RRPCs, uint8 pre_DEM_level, uint8 DEM_level,	uint8 NumOfIAparam, bool check_tile_array,bool Hemisphere,bool* tile_array,
			  CSize Limagesize,CSize Rimagesize,CSize LBRsize,CSize RBRsize,TransParam param,int total_count,float *ori_minmaxHeight,float *Boundary,int row_iter, int col_iter);
bool check_image_boundary(double **lrpc, double **rrpc, uint8 numofparam, float *rimageparam, F2DPOINT Lstartpos, F2DPOINT Rstartpos,
                          F2DPOINT pos_xy, float minH, float maxH, CSize Lsize, CSize Rsize, int H_template_size, int pyramid_step);
void RemoveFiles(char *save_path, char *lfilename, char *rfilename, int py_level, bool flag);
float MergeTiles(ProInfo info,int iter_row_end,int t_col_end, int buffer);


struct simplex *IMakeSimplex(int ndim1);
struct facelist *IMakeFacelist(int ndim);
struct nnweight *IMakeWeight();
int *IntVect();
void FreeVecti();
double *DoubleVect();
void FreeVectd();
int *IntMatrix();
void FreeMatrixi();
float *FloatMatrix();
void FreeMatrixf();
double *DoubleMatrix();
void FreeMatrixd();

void FindNebPts_F(NNXY *input, int row_size, int col_size, int grid, double minX, double minY, double maxX, double maxY, double X, double Y, int *numpts, int row_interval, int col_interval, int ndim1, char* path);
double *FindNebPts_F_M(NNXY *input, int row_size, int col_size, double grid, double minX, double minY, double maxX, double maxY, double X, double Y, int *numpts, int row_interval, int col_interval, int ndim1, char* path);
double FindNebPts_F_M_IDW(NNXY *input, int row_size, int col_size, double grid, double minX, double minY, double maxX, double maxY, double X, double Y, int *numpts, int row_interval, int col_interval, int ndim1, char* path);

double nnamain(int ndata, int ndim, int ndim1, double q1, double q2, int *check_find,char* path);
double nnamain_M(double *points,int ndata, int ndim, int ndim1, double q1, double q2, int *check_find, char* path);
double nnamain_M_R(double *points,int ndata, int ndim, int ndim1, double q1, double q2, int *check_find, char* path);

void NNA(char *save_path, char *datafile, int row, int col, int level, int iteration, float min_X, float min_Y, float grid_resolution, CSize Size_Grid2D);
void NNA_M(char *save_path, char *iterfile, int row_end, int col_end, float grid_resolution, float mt_grid_resolution, int buffer_clip, int Hemisphere);

void Envihdr_writer(char *filename, int col_size, int row_size, float grid_size, float minX, float maxY, int NS_flag, int data_type);
CSize Envihdr_reader(char *filename);
CSize Envihdr_reader_seedDEM(char *filename, float *minX, float *maxY, float *grid_size);
bool TFW_reader_seedDEM(char *filename, float *minX, float *maxY, float *grid_size);

//orthogeneration
void orthogeneration(char *ImageFilename, char *DEMFilename, char *Outputpath);
F2DPOINT OriginalToPyramid_single_ortho(F2DPOINT InCoord, F2DPOINT Startpos, uint8 Pyramid_step);
uint16 *Preprocessing_ortho(uint8 py_level, CSize *data_size, uint16 *subimg);
uint16* CreateImagePyramid_ortho(uint16* _input, CSize _img_size, int _filter_size, float _sigma);

void SetPySizes_ortho(CSize *data_size, CSize subsetsize, int level);

uint16 *subsetImage_ortho(TransParam transparam, double **RPCs, char *ImageFilename, 
						  float *subBoundary, float *minmaxHeight, F2DPOINT *startpos, char *subsetImage, CSize* subsetsize, bool *ret);

bool GetImageSize_ortho(char *filename, CSize *Imagesize);

bool GetsubareaImage_ortho(TransParam transparam, double **RPCs, char *ImageFilename, CSize *Imagesize,
						   float *subBoundary, float *minmaxHeight, int *cols, int *rows);

uint16 *Readtiff_ortho(char *filename, CSize Imagesize, int *cols, int *rows, CSize *data_size);
F2DPOINT GetObjectToImageRPC_single_ortho(double **_rpc, uint8 _numofparam, float *_imageparam, F3DPOINT _GP);
F2DPOINT* GetObjectToImageRPC_ortho(double **_rpc, uint8 _numofparam, float *_imageparam, uint16 _numofpts, F3DPOINT *_GP);

F2DPOINT *SetDEMGrid_ortho(float *Boundary, float Grid_x, float Grid_y, CSize Size_2D);

bool SetOrthoBoundary_ortho(CSize *Imagesize, float *Boundary, 
							double **RPCs, float gridspace, CSize DEM_size, float minX, float maxY, TransParam param);
float *LoadDEM_ortho(char *DEM_path, char* hdr_path);


double** OpenXMLFile_ortho(char* _filename, double* gsd_r, double* gsd_c);
CSize Envihdr_reader_ortho(char *filename);
void Envihdr_writer_ortho(char *filename, int col_size, int row_size, float grid_size, float minX, float maxY, int NS_flag, int data_type);
CSize Envihdr_reader_DEM_ortho(char *filename, float *minX, float *maxY, float *grid_size);
char* remove_ext_ortho(char* mystr);
F2DPOINT* wgs2ps_ortho(bool _bopenmp, TransParam _param, int _numofpts, F2DPOINT *_wgs);
F2DPOINT wgs2ps_single_ortho(TransParam _param, F2DPOINT _wgs);
F3DPOINT* wgs2ps_3D_ortho(bool _bopenmp, TransParam _param, int _numofpts, F3DPOINT *_wgs);
F2DPOINT* ps2wgs_ortho(bool _bopenmp, TransParam _param, int _numofpts, F2DPOINT *_ps);
F2DPOINT ps2wgs_single_ortho(TransParam _param, F2DPOINT _ps);
F3DPOINT* ps2wgs_3D_ortho(bool _bopenmp, TransParam _param, int _numofpts, F3DPOINT *_ps);
void SetTransParam_ortho(float minLat, bool *Hemisphere, TransParam *param);