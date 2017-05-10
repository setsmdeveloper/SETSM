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
	bool bHemisphere;
	// a and b are radius and eccentricity of WGS84 ellipsoid, repectively.
	double a, e;
	// phi_c and lambda_0 are the latitude of true scale of standard parallel and meridian along positive Y axis, respectively.
	double phi_c, lambda_0;
	double t_c, m_c;
	int pm;
	
	//UTM param
	double sa, sb, e2, e2cuadrada, c;
	int zone;
	char direction[10];
	int projection;
    int utm_zone;
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
	int roh_count;
	uint8 mag_tag;
	double INCC;
	double GNCC;
} NCCresult;

typedef struct UpdateGrid{
	double minHeight;
	double maxHeight;
	
	double t_minHeight;
	double t_maxHeight;
	
	double Height; //after blunder detection
	uint8 Matched_flag;
	double roh;
	uint8 anchor_flag;
	double Matched_height;//before blunder detection
	double ortho_ncc;
	double angle;
	uint16 false_h_count;
	double *false_h;
	
}UGRID;

typedef struct BlunderIP{
	uint8 Pyramid_step;
	CSize Size_Grid2D;
	double* Boundary;
	double gridspace;
	uint8 iteration;
	double Hinterval;
	bool height_check_flag;
}BL;

typedef struct ProjectInfo{
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
    
	uint8 SPnumber[2],NumOfTile_row, NumOfTile_col;	

	int threads_num;
	bool IsRA, IsSP, IsRR, IsSaveStep, Overall_DEM, Affine_RA, pre_DEMtif, check_tile_array;
} ProInfo;

typedef struct ArgumentInfo{
	int check_arg; // 0 : no input, 1: 3 input
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
	char Image1[500];
	char Image2[500];
	char Outputpath[500];
	char Outputpath_name[500];
	int Threads_num;
	double DEM_space;
	char seedDEMfilename[500];
	char metafilename[500];
	double seedDEMsigma;
	double minHeight;
	double maxHeight;
	int start_row;
	int end_row;
	int start_col;
	int end_col;
	double ra_line;
	double ra_sample;
	int RA_row;
	int RA_col;
	int tilesize;
	double Min_X, Max_X, Min_Y, Max_Y;
	
	int projection; //PS = 1, UTM = 2
    int utm_zone;
    int sensor_provider; //DG = 1, Pleiades = 2
    double image_resolution;
    int ortho_count;
    double overlap_length;
	
} ARGINFO;

#endif

#ifndef __VDEFS_H
#define __VDEFS_H

#ifndef NULL
#define NULL 0
#endif

#define DELETED -2

typedef struct tagFreenode
{
	struct tagFreenode * nextfree;
} Freenode ;


typedef struct tagFreelist
{
	Freenode * head;
	int nodesize;
} Freelist ;

typedef struct tagPoint
{
	double x ;
	double y ;
} Point ;

/* structure used both for sites and for vertices */

typedef struct tagSite
{
	Point coord ;
	int sitenbr ;
	int refcnt ;
} Site ;


typedef struct tagEdge
{
	double a, b, c ;
	Site * ep[2] ;
	Site * reg[2] ;
	int edgenbr ;
} Edge ;

#define le 0
#define re 1

typedef struct tagHalfedge
{
	struct tagHalfedge * ELleft ;
	struct tagHalfedge * ELright ;
	Edge * ELedge ;
	int ELrefcnt ;
	char ELpm ;
	Site * vertex ;
	double ystar ;
	struct tagHalfedge * PQnext ;
} Halfedge ;

/* edgelist.c */
void ELinitialize(void) ;
Halfedge * HEcreate(Edge *, int) ;
void ELinsert(Halfedge *, Halfedge *) ;
Halfedge * ELgethash(int) ;
Halfedge * ELleftbnd(Point *) ;
void ELdelete(Halfedge *) ;
Halfedge * ELright(Halfedge *) ;
Halfedge * ELleft(Halfedge *) ;
Site * leftreg(Halfedge *) ;
Site * rightreg(Halfedge *) ;
extern int ELhashsize ;
extern Site * bottomsite ;
extern Freelist hfl ;
extern Halfedge * ELleftend, * ELrightend, **ELhash ;

/* geometry.c */
void geominit(void) ;
Edge * bisect(Site *, Site *) ;
Site * intersect(Halfedge *, Halfedge *) ;
int right_of(Halfedge *, Point *) ;
void endpoint(Edge *, int, Site *) ;
double dist(Site *, Site *) ;
void makevertex(Site *) ;
void deref(Site *) ;
void ref(Site *) ;
extern double deltax, deltay ;
extern int nsites, nedges, sqrt_nsites, nvertices ;
extern Freelist sfl, efl ;

/* heap.c */
void PQinsert(Halfedge *, Site *, double) ;
void PQdelete(Halfedge *) ;
int PQbucket(Halfedge *) ;
int PQempty(void) ;
Point PQ_min(void) ;
Halfedge * PQextractmin(void) ;
void PQinitialize(void) ;
extern int PQmin, PQcount, PQhashsize ;
extern Halfedge * PQhash ;

/* getopt.c */
extern int getopt(int, char *const *, const char *);

/* memory.c */
void freeinit(Freelist *, int) ;
char *getfree(Freelist *) ;
void makefree(Freenode *, Freelist *) ;
char *myalloc(unsigned) ;
void free_all(void);

/* output.c */
void openpl(void) ;
void line(double, double, double, double) ;
void circle(double, double, double) ;
void range(double, double, double, double) ;
void out_bisector(Edge *) ;
void out_ep(Edge *) ;
void out_vertex(Site *) ;
void out_site(Site *) ;
void out_triple(Site *, Site *, Site *) ;
void plotinit(void) ;
void clip_line(Edge *) ;

/* voronoi.c */
void voronoi(Site *(*)(),UI3DPOINT* trilists) ;

#endif	


