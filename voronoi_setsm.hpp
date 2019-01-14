/* This software was adapted from the voronoi algorithm by Steven Fortune
(http://ect.bell-labs.com/who/sjf/) as modified by Derek Bradley
(http://zurich.disneyresearch.com/derekbradley/voronoi.html)

Reference: Steve J. Fortune (1987) A Sweepline Algorithm for Voronoi Diagrams,
Algorithmica 2, 153-174. */

#ifndef __VORONOI_H
#define __VORONOI_H

#include "Typedefine.hpp"
#include "tiff.h"
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

#ifdef __cplusplus
extern "C" {
#endif

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

/* heap.c */
void PQinsert(Halfedge *, Site *, double) ;
void PQdelete(Halfedge *) ;
int PQbucket(Halfedge *) ;
int PQempty(void) ;
Point PQ_min(void) ;
Halfedge * PQextractmin(void) ;
void PQinitialize(void) ;

/* main.c */
int scomp(const void * vs1, const void * vs2);
Site *nextone(void);
void readsites(D3DPOINT *ptslists,int numofpts);
Site *readone(void);
void initializeVoronoi(void);

/* memory.c */
void freeinit(Freelist *, int) ;
char *getfree(Freelist *) ;
void makefree(Freenode *, Freelist *) ;
char *myalloc(unsigned) ;
void free_all(void) ;

/* output.c */
void openpl(void) ;
void line(double, double, double, double) ;
void circle(double, double, double) ;
void range(double, double, double, double) ;
void out_bisector(Edge *) ;
void out_ep(Edge *) ;
void out_vertex(Site *) ;
void out_site(Site *) ;
void out_triple(Site *, Site *, Site *, int*) ;
void plotinit(void) ;
void clip_line(Edge *) ;

/* voronoi.c */
void voronoi(Site *(*)(void),UI3DPOINT*, int*) ;

#ifdef __cplusplus
}
#endif

#endif  


