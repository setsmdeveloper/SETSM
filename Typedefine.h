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

/*
 * Includes code derived from the voronoi algorithm by Steven Fortune
 * (http://ect.bell-labs.com/who/sjf/)
 * as modified by Derek Bradley
 * (http://zurich.disneyresearch.com/derekbradley/voronoi.html)
 *
 * Reference: Steve J. Fortune (1987) A Sweepline Algorithm for Voronoi Diagrams,
 * Algorithmica 2, 153-174.
 */

#include "tiff.h"

#ifndef _Typedefine_H_
#define _Typedefine_H_

#define PI 3.141592653589793
#define DegToRad PI/180
#define RadToDeg 180/PI

#ifndef bool
#define bool unsigned char
#define true 0x1
#define false 0x0
#endif

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
	double result0;
	double result1;
	double result2;
	double result3;
	double result4;
	
	double INCC;
	double GNCC;
	
	int roh_count;
	uint8 mag_tag;
	
} NCCresult;

typedef struct UpdateGrid{
	double minHeight;
	double maxHeight;
	
	double t_minHeight;
	double t_maxHeight;
	
	double Height; //after blunder detection
	double roh;
	double Matched_height;//before blunder detection
	double ortho_ncc;
	double angle;
//	double *false_h;

//	uint16 false_h_count;
	uint8 Matched_flag;
	uint8 anchor_flag;
}UGRID;

typedef struct BlunderIP{
	CSize Size_Grid2D;
	double gridspace;
	double Hinterval; 
	double* Boundary;
	uint8 Pyramid_step;
	uint8 iteration;
	bool height_check_flag;
}BL;

typedef struct ProjectInfo{
	double resolution;
	double DEM_resolution;
	double preDEM_space;
	double cal_boundary[4];
	double RA_param[2];
	double seedDEMsigma;
	
	double minHeight;
	double maxHeight;
	
	int start_row;
	int end_row;
	int start_col;
	int end_col;
	int threads_num;
	
	char LeftImagefilename[500];
	char RightImagefilename[500];
	char LeftRPCfilename[500];
	char RightRPCfilename[500];
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
    
	uint8 SPnumber[2],NumOfTile_row, NumOfTile_col;	
} ProInfo;

typedef struct ArgumentInfo{
	double DEM_space;
	double seedDEMsigma;
	double minHeight;
	double maxHeight;
	double ra_line;
	double ra_sample;
	double Min_X, Max_X, Min_Y, Max_Y;
	double image_resolution;
	double overlap_length;
	
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
    int sensor_provider; //DG = 1, Pleiades = 2
    int ortho_count;
	
    char Image1[500];
	char Image2[500];
	char Outputpath[500];
	char Outputpath_name[500];
	char seedDEMfilename[500];
	char metafilename[500];
	
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
    bool check_ortho;
    bool check_imageresolution;
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

#endif

