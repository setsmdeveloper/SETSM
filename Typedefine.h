/*
 * Copyright 2017 Myoung-Jong Noh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
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
} D3DPOINT;

typedef struct tagTransParam
{
	bool bHemisphere;
	// a and b are radius and eccentricity of WGS84 ellipsoid, repectively.
	float a, e;
	// phi_c and lambda_0 are the latitude of true scale of standard parallel and meridian along positive Y axis, respectively.
	float phi_c, lambda_0;
	double t_c, m_c;
	int pm;
    
    //UTM param
    double sa, sb, e2, e2cuadrada, c;
    int zone;
    char direction[10];
    int projection;
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
	float result0;
	float result1;
	float result2;
	float result3;
	float result4;
    int roh_count;
	uint8 mag_tag;
    float INCC;
    float GNCC;
} NCCresult;

typedef struct UpdateGrid{
	float minHeight;
	float maxHeight;
    
    float t_minHeight;
    float t_maxHeight;
    
	float Height; //after blunder detection
	uint8 Matched_flag;
	float roh;
	uint8 anchor_flag;
	float Matched_height;//before blunder detection
	float ortho_ncc;
    float angle;
    uint16 false_h_count;
    float *false_h;
    
}UGRID;

typedef struct BlunderIP{
	uint8 Pyramid_step;
	CSize Size_Grid2D;
	int* Boundary;
	float gridspace;
	uint8 iteration;
	float Hinterval;
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
	
	float resolution;
	float DEM_resolution;
	float preDEM_space;
	float cal_boundary[4];
	float RA_param[2];
    float seedDEMsigma;
	
	float minHeight;
	float maxHeight;
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
    
	char Image1[500];
	char Image2[500];
	char Outputpath[500];
    char Outputpath_name[500];
	int Threads_num;
	float DEM_space;
	char seedDEMfilename[500];
	char metafilename[500];
	float seedDEMsigma;
	float minHeight;
	float maxHeight;
	int start_row;
	int end_row;
	int start_col;
	int end_col;
	float ra_line;
	float ra_sample;
	int RA_row;
	int RA_col;
	int tilesize;
	int Min_X, Max_X, Min_Y, Max_Y;
	
    int projection;
    
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
    float x ;
    float y ;
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
    float a, b, c ;
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
    float ystar ;
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
float dist(Site *, Site *) ;
void makevertex(Site *) ;
void deref(Site *) ;
void ref(Site *) ;
extern float deltax, deltay ;
extern int nsites, nedges, sqrt_nsites, nvertices ;
extern Freelist sfl, efl ;

/* heap.c */
void PQinsert(Halfedge *, Site *, float) ;
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
void line(float, float, float, float) ;
void circle(float, float, float) ;
void range(float, float, float, float) ;
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


