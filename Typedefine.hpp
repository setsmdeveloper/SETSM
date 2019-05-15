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

#include "tiff.h"

#ifndef _Typedefine_H_
#define _Typedefine_H_

#define PI 3.141592653589793
#define DegToRad PI/180
#define RadToDeg 180/PI
#define UMToMM 0.001
#define MMToUM 1000
#define MaxImages 3

#ifndef bool
//#define bool unsigned char
#define true 0x1
#define false 0x0
#endif

enum SensorType {SB , AB};
enum SensorProvider {DG, PL};

typedef struct tagUI2DPoint
{
	uint32 m_X;
	uint32 m_Y;
} UI2DPOINT;

typedef struct tagUI3DPoint
{
	uint32 m_X;
	uint32 m_Y;
	uint32 m_Z;
} UI3DPOINT;


typedef struct tagF2DPoint
{
	float m_X;
	float m_Y;
} F2DPOINT;

typedef struct tagF3DPoint
{
	float m_X;
	float m_Y;
	float m_Z;
	uint8 flag;
} F3DPOINT;

typedef struct tagD2DPoint
{
	double m_X;
	double m_Y;
} D2DPOINT;

typedef struct tagD3DPoint
{
	double m_X;
	double m_Y;
	double m_Z;
    uint8 flag;
} D3DPOINT;

typedef struct tagTransParam
{
	double t_c, m_c;
	//UTM param
	double sa, sb, e2, e2cuadrada, c;
	// a and b are radius and eccentricity of WGS84 ellipsoid, repectively.
	double a, e;
	// phi_c and lambda_0 are the latitude of true scale of standard parallel and meridian along positive Y axis, respectively.
	double phi_c, lambda_0;
	
	int pm;
	int zone;
	int projection;
    int utm_zone;
	
	char direction[10];
	bool bHemisphere;
} TransParam;

typedef struct tagCSize
{
	unsigned int width;
	unsigned int height;
} CSize;

typedef struct tagNCCflag
{
	uint8 rotate_flag;
	uint8 multi_flag; 
	uint8 multi_flag_sum;
	uint8 inter_flag;
	uint8 weight_flag;
} NCCflag;

typedef struct tagNCCresult
{
	double result0; //first peak roh
	double result1; //second peak roh
	float result2; //first peak height
	float result3; //second peak height
    float GNCC;
	int result4; //peak count
    double max_WNCC;
    int max_WNCC_pos;
	
    int minHeight;
    int maxHeight;
    int NumOfHeight;
    bool check_height_change;
	//int roh_count;
} NCCresult;

typedef struct UpdateGrid{
	int minHeight;
	int maxHeight;
	
	float Height; //after blunder detection
	double roh;
	double ortho_ncc[MaxImages];
    double Mean_ortho_ncc;

    uint8 Matched_flag;
	uint8 anchor_flag;
    
}UGRID;

typedef struct tagVoxelinfo
{
    //float ANCC;
    float height;
    bool flag_cal;
    float INCC;
}VOXEL;

typedef struct LSFinfo{
    float lsf_std;
    unsigned char lsf_kernel;
}LSFINFO;

typedef struct BlunderIP{
	CSize Size_Grid2D;
	double gridspace;
	double Hinterval; 
	double* Boundary;
	uint8 Pyramid_step;
	uint8 iteration;
	bool height_check_flag;
}BL;

//Aerial frame Camera
typedef struct tagCameraInfo
{
    double m_focalLength;
    CSize m_ImageSize;
    double m_CCDSize;
    double m_ppx;
    double m_ppy;
} CAMERA_INFO;

typedef struct tagRotationMatrix
{
    double m11, m12, m13, m21, m22, m23, m31, m32, m33;
} RM;

typedef struct tagEO
{
    int strip_ID;
    int photo_ID;
    char path[500];
    double m_Xl;
    double m_Yl;
    double m_Zl;
    double m_Wl;
    double m_Pl;
    double m_Kl;
    RM m_Rm;
} EO;


typedef struct tagFrameInfo
{
    CAMERA_INFO m_Camera;
    int NumberofStip;
    int NumberofPhotos;
    int start_stripID;
    int end_stripID;
    EO* Photoinfo;
} FrameInfo;


typedef struct ProjectInfo{
	double resolution;
	double DEM_resolution;
	double preDEM_space;
	double cal_boundary[4];
	double RA_param[MaxImages][2];
	double seedDEMsigma;
	
	double minHeight;
	double maxHeight;
	double System_memory;
    
	int start_row;
	int end_row;
	int start_col;
	int end_col;
	int threads_num;
    int number_of_images;
    enum SensorType sensor_type; // 1 is for RFM (default), 2 is for Collinear Equation (Frame)
    uint8 pyramid_level;
    
	char Imagefilename[MaxImages][500];
	char RPCfilename[MaxImages][500];
	char save_filepath[500];
	char Outputpath_name[500];
	char tmpdir[500];
	char tile_info[500];
	char priori_DEM_tif[500];
	char metafilename[500];
	
	bool check_minH;
	bool check_maxH;
	bool check_tiles_SR;
	bool check_tiles_ER;
	bool check_tiles_SC;
	bool check_tiles_EC;
	bool check_gridonly;
	bool check_boundary;
	bool check_checktiff;
    bool check_ortho;
	bool IsRA, IsSP, IsRR, IsSaveStep, Overall_DEM, Affine_RA, pre_DEMtif, check_tile_array;
    bool check_Matchtag;
    bool check_selected_image[MaxImages];
    bool check_full_cal;
    
    //SGM test flag
    bool check_SNCC;
    bool check_updateheight;
    bool check_blunderdetection;
    bool check_NCCpeak;
    bool check_minTH;
    bool check_adaptive_P2;
    bool check_orthoblunder;
    bool check_8d;
    int SGM_py;
    
    FrameInfo frameinfo;
    
	uint8 SPnumber[2],NumOfTile_row, NumOfTile_col;	
} ProInfo;

typedef struct ArgumentInfo{
	double DEM_space;
	double seedDEMsigma;
	double minHeight;
	double maxHeight;
	double ra_line[MaxImages];
	double ra_sample[MaxImages];
	double Min_X, Max_X, Min_Y, Max_Y;
	double image_resolution;
	double overlap_length;
    double focal_length;
    double CCD_size;
    double System_memory;
    
	int check_arg; // 0 : no input, 1: 3 input
	int Threads_num;
	int start_row;
	int end_row;
	int start_col;
	int end_col;
	int RA_row;
	int RA_col;
	int tilesize;
	int projection; //PS = 1, UTM = 2
    int utm_zone;
    enum SensorType sensor_type; // SB is for RFM (default), AB is for Collinear Equation (Frame)
    enum SensorProvider sensor_provider; //DG = DG, Pleiades = PL if sensor_type = 1
    int ortho_count;
    int RA_only;
    int number_of_images; // 2 is for stereo (default), n is for multi more than 3
    uint8 pyramid_level;
    
    char Image[MaxImages][500];
	char Outputpath[500];
	char Outputpath_name[500];
	char seedDEMfilename[500];
	char metafilename[500];
    char EO_Path[500];
	
	bool check_DEM_space;
	bool check_Threads_num;
	bool check_seeddem;
	bool check_minH;
	bool check_maxH;
	bool check_tiles_SR;
	bool check_tiles_ER;
	bool check_tiles_SC;
	bool check_tiles_EC;
	bool check_RA_line;
	bool check_RA_sample;
	bool check_gridonly; 
	bool check_RA_tileR;
	bool check_RA_tileC;
	bool check_tilesize;
	bool check_boundary;
	bool check_checktiff;
	bool check_RA_only;
    bool check_ortho;
    bool check_imageresolution;
    bool check_LSF;
    bool check_LSF_DEM;
    bool check_LSFDEMpath;
    int check_LSF2;
    bool check_Matchtag;
    bool check_EO;
    bool check_fl;
    bool check_ccd;
    bool check_full_cal;
    
    //SGM test flag
    bool check_SNCC;
    bool check_updateheight;
    bool check_blunderdetection;
    bool check_NCCpeak;
    bool check_minTH;
    bool check_adaptive_P2;
    bool check_orthoblunder;
    bool check_8d;
    int SGM_py;
    
} ARGINFO;

typedef struct tagImageInfo
{
    float Mean_sun_azimuth_angle;
    float Mean_sun_elevation;
    float Mean_sat_azimuth_angle;
    float Mean_sat_elevation;
    float Intrack_angle;
    float Crosstrack_angle;
    float Offnadir_angle;
    float cloud;
    float Image_ori;
    float Image_ori_azi;
    float dx,dy,f;
    float UL[3], UR[3],LR[3],LL[3];
    float convergence_angle;
	
	int month;
    int date;
    int year;
    int scandirection;
    
    char filename[500];
    char Ofilename[500];
    char imagetime[500];
    char SatID[500];
} ImageInfo;

typedef struct tagImageGSD
{
    float row_GSD;
    float col_GSD;
    float pro_GSD;
} ImageGSD;

typedef struct tagBandInfo
{
    float abscalfactor;
    float effbw;
    float tdi;
} BandInfo;

typedef struct
{
    long int ncols;
    long int nrows;
    long double **val;
    long double *data;
} GMA_double;

#endif

