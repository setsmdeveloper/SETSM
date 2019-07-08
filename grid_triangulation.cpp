/***************************************
 * The triangulation algorithm in the code was adapted from the methods in these 
 * two papers:
 * [1] Wenzhou Wu, Yikang Rui, Fenzhen Su, Liang Cheng & Jiechen Wang (2014) 
 * Novel parallel algorithm for constructing Delaunay triangulation based on a
 * twofold-divide-and-conquer scheme, GIScience & Remote Sensing, 51:5, 537-554, DOI:
 * 10.1080/15481603.2014.946666
 * [2] Leonidas Guibas & Jorge Stolfi (1985) Primitives for the Manipulation of 
 * General Subdivisions and the Computation of Voronoi Diagrams, ACM Transactions on 
 * Graphics, 51:2, 74-123.
 ***************************************/

#include <iostream>
#include <algorithm>
#include <vector>
#include <omp.h>

#include "grid_triangulation.hpp"

//////////////////////////////////////
// HELPERS ///////////////////////////
//////////////////////////////////////

/* Positive if a, b, c counter-clockwise
 * Negative if a, b, c clockwise
 * Zero if a, b, c collinear
 */
inline INT64 Orientation(const GridPoint &a, const GridPoint &b, const GridPoint &c)
{
	INT64 d_11 = a.col - c.col;
	INT64 d_21 = b.col - c.col;
	INT64 d_12 = a.row - c.row;
	INT64 d_22 = b.row - c.row;
	
	return d_11 * d_22 - d_12 *d_21;
}

/* Assumes a, b, c in counter-clockwise order
 * Then	positive if d in circle,
 * 	negative if d outside circle,
 * 	zero if d on circle
 */
inline INT128 InCircle(const GridPoint &a, const GridPoint &b, const GridPoint &c, const GridPoint &d)
{
	INT64 d_11 = a.col - d.col;
	INT64 d_12 = a.row - d.row;
	INT64 d_13 = d_11 * d_11 + d_12 * d_12;
	INT64 d_21 = b.col - d.col;
	INT64 d_22 = b.row - d.row;
	INT64 d_23 = d_21 * d_21 + d_22 * d_22;
	INT64 d_31 = c.col - d.col;
	INT64 d_32 = c.row - d.row;
	INT64 d_33 = d_31 * d_31 + d_32 * d_32;

	INT128 bc = d_21 * d_32 - d_31 * d_22;
	INT128 ac = d_11 * d_32 - d_31 * d_12;
	INT128 ab = d_11 * d_22 - d_21 * d_12;

	return d_13 * bc - d_23 * ac + d_33 * ab;
}

inline bool LessThanPtrXY(GridPoint *a, GridPoint *b)
{
	if (a && b) return LessThanXY(*a, *b);
	else return false;
}

inline bool LessThanPtrYX(GridPoint *a, GridPoint *b)
{
	if (a && b) return LessThanYX(*a, *b);
	else return false;
}

//////////////////////////////////////
// GridTriangulation Implementation //
//////////////////////////////////////

template <typename GridType, typename IterType>
GridTriangulation<GridType, IterType>::GridTriangulation(INDEX width, INDEX height)
{
	this->width = width;
	this->height = height;
	this->grid = 0;
	this->edge_list = 0;
}

template <typename GridType, typename IterType>
GridTriangulation<GridType, IterType>::~GridTriangulation()
{
	delete this->grid;
	delete this->edge_list;
}


template <typename GridType, typename IterType>
Grid<GridType, IterType> *GridTriangulation<GridType, IterType>::getGrid()
{
	return this->grid;
}


template <typename GridType, typename IterType>
Edge *GridTriangulation<GridType, IterType>::AddEdgeAndTwin(const GridPoint &orig, const GridPoint &dest)
{
	// Get edges from edge list
	Edge *e = this->edge_list->GetNewEdge();
	Edge *et = this->edge_list->GetNewEdge(); 

	// Set next/prev appropriately
	e->twin = (e->oprev = (e->dnext = et));
	et->twin = (et->oprev = (et->dnext = e));

	// Tie to orig/dest
	e->orig = orig;
	et->orig = dest;
	
	return e;
}

template <typename GridType, typename IterType>
void GridTriangulation<GridType, IterType>::RemoveEdgeAndTwin(Edge &edge)
{
	// Get edge, its twin, and orig/dest
	Edge *e = &edge;
	Edge *et = e->twin;	
	GridPoint &orig = e->orig;
	GridPoint &dest = et->orig;

	// Update next/prev edges
	e->oprev->dnext = et->dnext;
	e->dnext->oprev = et->oprev;
	et->oprev->dnext = e->dnext;
	et->dnext->oprev = e->oprev;

	// Free edge and twin from edge_list
	// so memory can be used for new edges
	this->edge_list->RemoveEdge(*e);
	this->edge_list->RemoveEdge(*et);
}

template <typename GridType, typename IterType>
void GridTriangulation<GridType, IterType>::Weld(Edge &in, Edge &out)
{
	// Set edges before and after
	// in (counterclockwise)	
	Edge *prev = out.oprev;
	Edge *curr = &in;
	Edge *next = &out;

	// Update next/prev edges
	next->oprev = curr;
	curr->dnext = next;
	curr->twin->oprev = prev;
	prev->dnext = curr->twin;
}
	
template <typename GridType, typename IterType>
Edge *GridTriangulation<GridType, IterType>::Bridge(Edge &start, Edge &end)
{
	Edge *in = &start;
	Edge *out = &end;
	GridPoint &orig = in->twin->orig;
	GridPoint &dest = out->orig;

	// Create bridging edge and its twin
	Edge *e = this->AddEdgeAndTwin(orig, dest);
	Edge *et = e->twin;

	// Make appropriate welds to connect the
	// new edges to the rest of the graph
	this->Weld(*et, *(in->dnext));
	this->Weld(*e, *out);

	return e;
}

template <typename GridType, typename IterType>
void GridTriangulation<GridType, IterType>::Triangulate(GridPoint *points[], size_t num_points)
{
	// Initialize memory for edges by creating an EdgeList
	// and run the triangulation routine
	// Note: ex contains some of the edges of the convex hull,
	// not usefull for result of simple triangulation, so it is
	// discarded
	#pragma omp parallel
	{
		#pragma omp single
		{
			//double b = omp_get_wtime();
			this->edge_list = new EdgeList(num_points);
			ExtremeEdges *ex = this->TriangulateHorizontalThreaded(points, num_points);
			delete ex;
			//double e = omp_get_wtime();
			//fprintf(stderr, "\tTriangulation took %f\n", e - b);
		}
	}

	//double b = omp_get_wtime();
	this->grid = new Grid<GridType, IterType>(this->width, this->height, points, num_points);
	this->edge_list->SetGrid(this->grid);
	//double e = omp_get_wtime();
	//fprintf(stderr, "\tgrid setting took %f\n", e - b);
}

template <typename GridType, typename IterType>
ExtremeEdges *GridTriangulation<GridType, IterType>::TriangulateHorizontalThreaded(GridPoint *points[], size_t num_points)
{
	// Base Cases /////////////////////////////////
	if (num_points < THREADED_CUTOFF)
	{
		return this->TriangulateHorizontalSerial(points, num_points);
	}
	else if (num_points < 2)
	{
		return 0;
	}
	else if (num_points == 2)
	{
		return this->Triangulate2(points);
	}
	else if (num_points == 3)
	{
		return this->Triangulate3(points);
	}
	
	// Recursive Case /////////////////////////////
	
	// Partition lexicographically (compare x/col, then y/row)
	// on median of the list
	size_t median = num_points / 2;
	std::nth_element(points, points + median, points + num_points, LessThanPtrXY);	

	// Perform recursive triangulations of list halves
	// Note alternating direction of cuts for better expected runtime
	GridTriangulation *g_left = new GridTriangulation(this->width, this->height);
	GridTriangulation *g_right = new GridTriangulation(this->width, this->height);
	this->edge_list->SplitEdgeList(g_left->edge_list, g_right->edge_list, median);

	ExtremeEdges *left_ex;
	ExtremeEdges *right_ex;
	#pragma omp task shared(left_ex)
	{
		left_ex = g_left->TriangulateVerticalThreaded(points, median);
	}
	#pragma omp task shared(right_ex)
	{
		right_ex = g_right->TriangulateVerticalThreaded(points + median, num_points - median);
	}
	#pragma omp taskwait

	this->edge_list->MergeEdgeLists(*(g_left->edge_list), *(g_right->edge_list));
	delete g_left;
	delete g_right;

	ExtremeEdges *ex = new ExtremeEdges();

	// Use extreme edges from recursive triangulations to determine
	// the lower common tangent (lct) of combined convex hull
	Edge *left_edge = left_ex->right_edge_cw;
	Edge *right_edge = right_ex->left_edge_ccw;
	while (1)
	{
		if (Orientation(left_edge->orig, left_edge->twin->orig, right_edge->orig) > 0) left_edge = left_edge->twin->oprev->twin;
		else if (Orientation(right_edge->orig, right_edge->twin->orig, left_edge->orig) < 0) right_edge = right_edge->dnext;
		else break;
	}
	Edge *lct = this->Bridge(*(right_edge->oprev), *(left_edge->twin->dnext));

	// Check if lct 'covers' extreme edges from either sub-triangulation	
	if (left_edge->orig == left_ex->left_edge_ccw->orig) left_ex->left_edge_ccw = lct->twin;
	if (right_edge->orig == right_ex->right_edge_cw->orig) right_ex->right_edge_cw = lct;

	// Perform merge of two sub-triangulations, ending with the
	// creation of the upper common tangent (uct) of the combined
	// convex hull
	Edge *new_base = lct;
	Edge *uct = 0;
	while (new_base != 0)
	{
		uct = new_base;
		new_base = this->NextCrossEdge(uct);
	}

	// Set the left/right extreme edges of the combined
	// convex hull using the results from the sub-triangulations
	ex->left_edge_ccw = left_ex->left_edge_ccw;
	ex->right_edge_cw = right_ex->right_edge_cw;
	
	// Set the bottom extreme edge of the combined
	// convex hull using the lct as a starting guess
	Edge *temp = lct->twin;
	while (LessThanYX(temp->orig, temp->twin->orig)) temp = temp->oprev;
	temp = temp->dnext;
	while (LessThanYX(temp->twin->orig, temp->orig)) temp = temp->dnext;
	ex->bottom_edge_ccw = temp;

	// Set the top extreme edge of the combined
	// convex hull using the uct as a starting guess
	temp = uct->twin;
	while (LessThanYX(temp->twin->orig, temp->orig)) temp = temp->twin->dnext->twin;
	temp = temp->twin->oprev->twin;
	while (LessThanYX(temp->orig, temp->twin->orig)) temp = temp->twin->oprev->twin;
	ex->top_edge_cw = temp;

	delete left_ex;
	delete right_ex;
	return ex;
}

template <typename GridType, typename IterType>
ExtremeEdges *GridTriangulation<GridType, IterType>::TriangulateVerticalThreaded(GridPoint *points[], size_t num_points)
{
	// Base Cases /////////////////////////////////
	if (num_points < THREADED_CUTOFF)
	{
		return this->TriangulateVerticalSerial(points, num_points);
	}
	else if (num_points < 2)
	{
		return 0;
	}
	else if (num_points == 2)
	{
		return this->Triangulate2(points);
	}
	else if (num_points == 3)
	{
		return this->Triangulate3(points);
	}

	// Recursive Case /////////////////////////////

	// Partition (compare y/row, then -x/col) points on
	// median of the list
	size_t median = num_points / 2;
	std::nth_element(points, points + median, points + num_points, LessThanPtrYX);

	// Perform recursive triangulations of list halves
	// Note alternating direction of cuts for better expected runtime
	GridTriangulation *g_bottom = new GridTriangulation(this->width, this->height);
	GridTriangulation *g_top = new GridTriangulation(this->width, this->height);
	this->edge_list->SplitEdgeList(g_bottom->edge_list, g_top->edge_list, median);

	ExtremeEdges *bottom_ex;
	ExtremeEdges *top_ex;
	#pragma omp task shared(bottom_ex)
	{
		bottom_ex = g_bottom->TriangulateHorizontalThreaded(points, median);
	}
	#pragma omp task shared(top_ex)
	{
		top_ex = g_top->TriangulateHorizontalThreaded(points + median, num_points - median);
	}
	#pragma omp taskwait

	this->edge_list->MergeEdgeLists(*(g_bottom->edge_list), *(g_top->edge_list));
	delete g_bottom;
	delete g_top;

	ExtremeEdges *ex = new ExtremeEdges();

	// Use extreme edges from recursive triangulations to determine
	// the right common tangent (rct) of combined convex hull
	Edge *bottom_edge = bottom_ex->top_edge_cw;
	Edge *top_edge = top_ex->bottom_edge_ccw;
	while (1)
	{
		if (Orientation(bottom_edge->orig, bottom_edge->twin->orig, top_edge->orig) > 0) bottom_edge = bottom_edge->twin->oprev->twin;
		else if (Orientation(top_edge->orig, top_edge->twin->orig, bottom_edge->orig) < 0) top_edge = top_edge->dnext;
		else break;
	}
	Edge *rct = this->Bridge(*(top_edge->oprev), *(bottom_edge->twin->dnext));
	
	// Check if rct 'covers' extreme edges from either sub-triangulation	
	if (bottom_edge->orig == bottom_ex->bottom_edge_ccw->orig) bottom_ex->bottom_edge_ccw = rct->twin;
	if (top_edge->orig == top_ex->top_edge_cw->orig) top_ex->top_edge_cw = rct;

	// Perform merge of two sub-triangulations, ending with the
	// creation of the left common tangent (lct) of the combined
	// convex hull
	Edge *new_base = rct;
	Edge *lct = 0;
	while (new_base != 0)
	{
		lct = new_base;
		new_base = this->NextCrossEdge(lct);
	}

	// Set the bottom/top extreme edges of the combined
	// convex hull using the results from the sub-triangulations
	ex->bottom_edge_ccw = bottom_ex->bottom_edge_ccw;
	ex->top_edge_cw = top_ex->top_edge_cw;

	// Set the right extreme edge of the combined
	// convex hull using the rct as a starting guess
	Edge *temp = rct;
	while (LessThanXY(temp->twin->orig, temp->orig)) temp = temp->twin->dnext->twin;
	temp = temp->twin->oprev->twin;
	while (LessThanXY(temp->orig, temp->twin->orig)) temp = temp->twin->oprev->twin;
	ex->right_edge_cw = temp;

	// Set the left extreme edge of the combined
	// convex hull using the lct as a starting guess
	temp = lct;
	while (LessThanXY(temp->orig, temp->twin->orig)) temp = temp->oprev;
	temp = temp->dnext;
	while (LessThanXY(temp->twin->orig, temp->orig)) temp = temp->dnext;
	ex->left_edge_ccw = temp;

	delete bottom_ex;
	delete top_ex;
	return ex;
}

template <typename GridType, typename IterType>
ExtremeEdges *GridTriangulation<GridType, IterType>::TriangulateHorizontalSerial(GridPoint *points[], size_t num_points)
{
	// Base Cases /////////////////////////////////
	if (num_points < 2)
	{
		return 0;
	}
	else if (num_points == 2)
	{
		return this->Triangulate2(points);
	}
	else if (num_points == 3)
	{
		return this->Triangulate3(points);
	}
	
	// Recursive Case /////////////////////////////
	
	// Partition lexicographically (compare x/col, then y/row)
	// on median of the list
	size_t median = num_points / 2;
	std::nth_element(points, points + median, points + num_points, LessThanPtrXY);	

	// Perform recursive triangulations of list halves
	// Note alternating direction of cuts for better expected runtime
	ExtremeEdges *left_ex = this->TriangulateVerticalSerial(points, median);
	ExtremeEdges *right_ex = this->TriangulateVerticalSerial(points + median, num_points - median);
	ExtremeEdges *ex = new ExtremeEdges();

	// Use extreme edges from recursive triangulations to determine
	// the lower common tangent (lct) of combined convex hull
	Edge *left_edge = left_ex->right_edge_cw;
	Edge *right_edge = right_ex->left_edge_ccw;
	while (1)
	{
		if (Orientation(left_edge->orig, left_edge->twin->orig, right_edge->orig) > 0) left_edge = left_edge->twin->oprev->twin;
		else if (Orientation(right_edge->orig, right_edge->twin->orig, left_edge->orig) < 0) right_edge = right_edge->dnext;
		else break;
	}
	Edge *lct = this->Bridge(*(right_edge->oprev), *(left_edge->twin->dnext));

	// Check if lct 'covers' extreme edges from either sub-triangulation	
	if (left_edge->orig == left_ex->left_edge_ccw->orig) left_ex->left_edge_ccw = lct->twin;
	if (right_edge->orig == right_ex->right_edge_cw->orig) right_ex->right_edge_cw = lct;

	// Perform merge of two sub-triangulations, ending with the
	// creation of the upper common tangent (uct) of the combined
	// convex hull
	Edge *new_base = lct;
	Edge *uct = 0;
	while (new_base != 0)
	{
		uct = new_base;
		new_base = this->NextCrossEdge(uct);
	}

	// Set the left/right extreme edges of the combined
	// convex hull using the results from the sub-triangulations
	ex->left_edge_ccw = left_ex->left_edge_ccw;
	ex->right_edge_cw = right_ex->right_edge_cw;
	
	// Set the bottom extreme edge of the combined
	// convex hull using the lct as a starting guess
	Edge *temp = lct->twin;
	while (LessThanYX(temp->orig, temp->twin->orig)) temp = temp->oprev;
	temp = temp->dnext;
	while (LessThanYX(temp->twin->orig, temp->orig)) temp = temp->dnext;
	ex->bottom_edge_ccw = temp;

	// Set the top extreme edge of the combined
	// convex hull using the uct as a starting guess
	temp = uct->twin;
	while (LessThanYX(temp->twin->orig, temp->orig)) temp = temp->twin->dnext->twin;
	temp = temp->twin->oprev->twin;
	while (LessThanYX(temp->orig, temp->twin->orig)) temp = temp->twin->oprev->twin;
	ex->top_edge_cw = temp;

	delete left_ex;
	delete right_ex;
	return ex;
}

template <typename GridType, typename IterType>
ExtremeEdges *GridTriangulation<GridType, IterType>::TriangulateVerticalSerial(GridPoint *points[], size_t num_points)
{
	// Base Cases /////////////////////////////////
	if (num_points < 2)
	{
		return 0;
	}
	else if (num_points == 2)
	{
		return this->Triangulate2(points);
	}
	else if (num_points == 3)
	{
		return this->Triangulate3(points);
	}

	// Recursive Case /////////////////////////////

	// Partition (compare y/row, then -x/col) points on
	// median of the list
	size_t median = num_points / 2;
	std::nth_element(points, points + median, points + num_points, LessThanPtrYX);

	// Perform recursive triangulations of list halves
	// Note alternating direction of cuts for better expected runtime
	ExtremeEdges *bottom_ex = this->TriangulateHorizontalSerial(points, median);
	ExtremeEdges *top_ex = this->TriangulateHorizontalSerial(points + median, num_points - median);
	ExtremeEdges *ex = new ExtremeEdges();

	// Use extreme edges from recursive triangulations to determine
	// the right common tangent (rct) of combined convex hull
	Edge *bottom_edge = bottom_ex->top_edge_cw;
	Edge *top_edge = top_ex->bottom_edge_ccw;
	while (1)
	{
		if (Orientation(bottom_edge->orig, bottom_edge->twin->orig, top_edge->orig) > 0) bottom_edge = bottom_edge->twin->oprev->twin;
		else if (Orientation(top_edge->orig, top_edge->twin->orig, bottom_edge->orig) < 0) top_edge = top_edge->dnext;
		else break;
	}
	Edge *rct = this->Bridge(*(top_edge->oprev), *(bottom_edge->twin->dnext));
	
	// Check if rct 'covers' extreme edges from either sub-triangulation	
	if (bottom_edge->orig == bottom_ex->bottom_edge_ccw->orig) bottom_ex->bottom_edge_ccw = rct->twin;
	if (top_edge->orig == top_ex->top_edge_cw->orig) top_ex->top_edge_cw = rct;

	// Perform merge of two sub-triangulations, ending with the
	// creation of the left common tangent (lct) of the combined
	// convex hull
	Edge *new_base = rct;
	Edge *lct = 0;
	while (new_base != 0)
	{
		lct = new_base;
		new_base = this->NextCrossEdge(lct);
	}

	// Set the bottom/top extreme edges of the combined
	// convex hull using the results from the sub-triangulations
	ex->bottom_edge_ccw = bottom_ex->bottom_edge_ccw;
	ex->top_edge_cw = top_ex->top_edge_cw;

	// Set the right extreme edge of the combined
	// convex hull using the rct as a starting guess
	Edge *temp = rct;
	while (LessThanXY(temp->twin->orig, temp->orig)) temp = temp->twin->dnext->twin;
	temp = temp->twin->oprev->twin;
	while (LessThanXY(temp->orig, temp->twin->orig)) temp = temp->twin->oprev->twin;
	ex->right_edge_cw = temp;

	// Set the left extreme edge of the combined
	// convex hull using the lct as a starting guess
	temp = lct;
	while (LessThanXY(temp->orig, temp->twin->orig)) temp = temp->oprev;
	temp = temp->dnext;
	while (LessThanXY(temp->twin->orig, temp->orig)) temp = temp->dnext;
	ex->left_edge_ccw = temp;

	delete bottom_ex;
	delete top_ex;
	return ex;
}

template <typename GridType, typename IterType>
ExtremeEdges *GridTriangulation<GridType, IterType>::Triangulate2(GridPoint *points[])
{
	ExtremeEdges *ex = new ExtremeEdges();
	
	// Create edge and twin between two points
	GridPoint a = *(points[0]);
	GridPoint b = *(points[1]);
	Edge *e = this->AddEdgeAndTwin(a, b);

	// Check lexicographic order (compare x/col, then y/row)
	// of two points to determine left/right extreme edges
	if (LessThanXY(a, b))
	{
		ex->left_edge_ccw = e;
		ex->right_edge_cw = e->twin;
	}
	else
	{
		ex->left_edge_ccw = e->twin;
		ex->right_edge_cw = e;
	}

	// Check lexicographic order (compare y/row, then -x/col)
	// of two points to determine top/bottom extreme edges
	if (LessThanYX(a, b))
	{
		ex->bottom_edge_ccw = e;
		ex->top_edge_cw = e->twin;
	}
	else
	{
		ex->bottom_edge_ccw = e->twin;
		ex->top_edge_cw = e;
	}

	return ex;
}

template <typename GridType, typename IterType>
ExtremeEdges *GridTriangulation<GridType, IterType>::Triangulate3(GridPoint *points[])
{
	ExtremeEdges *ex = new ExtremeEdges();

	// Sort points lexicographically (compare x/col, then y/row)
	//std::nth_element(points, points + 1, points + 3, LessThanPtrXY);
	std::sort(points, points + 3, LessThanPtrXY);
	GridPoint a = *(points[0]);
	GridPoint b = *(points[1]);
	GridPoint c = *(points[2]);

	// Add edges/twins between a,b and b,c
	Edge *e1 = this->AddEdgeAndTwin(a, b);
	Edge *e2 = this->AddEdgeAndTwin(b, c);
	this->Weld(*e1, *e2);

	// Check sign of area of triangle a,b,c
	// to:
	// 	a) see if third edge is necessary
	// 	b) determine right/left extreme edges 
	Edge *e3;
	INT64 signed_area = Orientation(a, b, c);
	if (signed_area > 0) // a, b, c is a counterclockwise-oriented triangle
	{
		e3 = this->Bridge(*e2, *e1);
		ex->left_edge_ccw = e1;
		ex->right_edge_cw = e2->twin;
	}
	else if (signed_area < 0) // a, b, c is a clockwise-oriented triangle
	{
		e3 = this->Bridge(*e2, *e1);
		ex->left_edge_ccw = e3->twin;
		ex->right_edge_cw = e3;
	}
	else // a, b, c are collinear
	{
		ex->left_edge_ccw = e1;
		ex->right_edge_cw = e2->twin;
	}

	// Sort points (compare y/row, -x/col)
	//std::nth_element(points, points + 1, points + 3, LessThanPtrYX);
	std::sort(points, points + 3, LessThanPtrYX);
	a = *(points[0]);
	b = *(points[1]);
	c = *(points[2]);

	// Collect edges previously created so that
	// e1 is the edge from a to b
	// and e2, e3 follow in order

	//Edge *e = this->GetEdgeOut(a);
	
	Edge *e = (e1->orig == a) ? e1 : (e2->orig == a ? e2 : e2->twin);
	e1 = (e->twin->orig == b) ? e : e->twin->dnext;
	e2 = e1->dnext;
	e3 = e2->dnext;

	// Check sign of area of triangle a,b,c
	// to determine top/bottom extreme edges
	signed_area = Orientation(a, b, c);
	if (signed_area >= 0)
	{
		ex->bottom_edge_ccw = e1;
		ex->top_edge_cw = e2->twin;
	}
	else
	{
		ex->bottom_edge_ccw = e3->twin;
		ex->top_edge_cw = e3;
	}

	return ex;
}

template <typename GridType, typename IterType>
Edge *GridTriangulation<GridType, IterType>::NextCrossEdge(Edge *base)
{
	// Determine first valid candidate for next left-to-right
	// bridge in merge operation, setting valid_l to false
	// if no such candidates exist
	Edge *l_cand = base->dnext;
	bool valid_l = Orientation(base->orig, l_cand->orig, l_cand->twin->orig) < 0;
	if (valid_l)
	{
		Edge *next_cand = l_cand->twin->dnext;
		while (InCircle(l_cand->orig, base->orig, l_cand->twin->orig, next_cand->twin->orig) > 0)
		{
			this->RemoveEdgeAndTwin(*l_cand);
			l_cand = next_cand;
			next_cand = l_cand->twin->dnext;
		}
	}

	// Determine first valid candidate for next right-to-left
	// bridge in merge operation, setting valid_r to false
	// if no such candidates exist
	Edge *r_cand = base->oprev->twin;
	bool valid_r = Orientation(base->twin->orig, base->orig, r_cand->twin->orig) > 0;
	if (valid_r)
	{
		Edge *next_cand = r_cand->oprev->twin;
		while (InCircle(base->twin->orig, base->orig, r_cand->twin->orig, next_cand->twin->orig) > 0)
		{
			this->RemoveEdgeAndTwin(*r_cand);
			r_cand = next_cand;
			next_cand = r_cand->oprev->twin;
		}
	}

	// Return null if no candidates were found,
	// otherwise pick valid candidate.
	// Existence of delaunay triangulation guarentees
	// candidate, and in for general position points,
	// it is unique. In special case where both
	// left-to-right and right-to-left candidates
	// are valid, the code arbitrarily select the
	// left-to-right candidate
	if (!valid_l && !valid_r) return 0;
	else if (!valid_l || (valid_r && (InCircle(l_cand->twin->orig, l_cand->orig, r_cand->orig, r_cand->twin->orig) > 0))) return (this->Bridge(*base, *(r_cand->twin)))->twin;
	else return (this->Bridge(*l_cand, *base))->twin;
}

/* Assumes point is included in full
 * triangulation.
 * Returns 0 if point is not on
 * convex hull, returns the clockwise
 * edge out of p on the convex hull
 * otherwise.
 */
template <typename GridType, typename IterType>
Edge *GridTriangulation<GridType, IterType>::OnConvexHull(const GridPoint &p)
{
	// Check that point is actually part of triangulation,
	// that is, has an edge out
	Edge *e = this->GetEdgeOut(p);
	if (e == 0) return 0;

	// Iterate over adjacent edges out of p
	// and check for reflex angles
	Edge *f = e;
	do
	{
		GridPoint &a = f->twin->orig;
		GridPoint &b = f->twin->dnext->twin->orig;
		GridPoint &c = f->orig;

		if (Orientation(a, b, c) <= 0) return f;
	} while ((f = f->twin->dnext) != e);

	return 0;
}

template <typename GridType, typename IterType>
void GridTriangulation<GridType, IterType>::TriangulateEmptyPolygon(Edge &edge)
{
	Edge *e = &edge;
	int n = 1;
	while ((e = e->dnext) != &edge) n++;

	// Base Case
	if (n <= 3) return;

	// Recursive Case
	std::vector<Edge *> edges;
	edges.reserve(n);
	for (int i = 0; i < n; ++i)
	{
		edges.push_back(e);
		e = e->dnext;
	}

	// Check if tri p_1, p_0, p_i is delaunay,
	// for i = 2, ..., n-1
	for (int i = 2; i < n; ++i)
	{
		bool is_delaunay = true;
		for (int j = 2; j < n; ++j)
		{
			if ((j != i) && (InCircle(e->twin->orig, e->orig, edges[i]->orig, edges[j]->orig) > 0))
			{
				is_delaunay = 0;
				break;
			}
		}

		if (is_delaunay)
		{
			if (i == 2)
			{
				e = this->Bridge(*(edges[1]), *(edges[0]));
				this->TriangulateEmptyPolygon(*(e->twin));
			}
			else if (i == n - 1)
			{
				e = this->Bridge(*(edges[0]), *(edges[n-1]));
				this->TriangulateEmptyPolygon(*(e->twin));
			}
			else
			{
				Edge *e1 = this->Bridge(*(edges[i-1]), *(edges[1]));
				Edge *e2 = this->Bridge(*(edges[n-1]), *(edges[i]));
				this->TriangulateEmptyPolygon(*e1);
				this->TriangulateEmptyPolygon(*e2);
			}
			break;
		}
	}
}

template <typename GridType, typename IterType>
void GridTriangulation<GridType, IterType>::TriangulateBorder(Edge& e, const GridPoint& p)
{
	for (Edge *e1 = &e; e1->orig != p; e1 = e1->dnext)
	{
		GridPoint a = e1->orig;
		for (Edge *e2 = e1->dnext; e2->orig != p; e2 = e2->dnext)
		{
			GridPoint b = e2->twin->orig;
			bool on_hull = true;
			for (Edge *e3 = e1; e3 != e2; e3 = e3->dnext) on_hull &= (Orientation(a, b, e3->twin->orig) > 0);
			for (Edge *e3 = e2; e3->orig != p; e3 = e3->dnext) on_hull &= (Orientation(a, b, e3->twin->orig) >= 0);
			if (on_hull)
			{
				this->TriangulateBorder(*(e2->dnext), p);
				this->TriangulateEmptyPolygon(*(this->Bridge(*e2, *e1)));
				return;
			}
		}
	}
}

template <typename GridType, typename IterType>
void GridTriangulation<GridType, IterType>::RemovePointAndRetriangulate(const GridPoint &p)
{
	Edge *e_p = this->GetEdgeOut(p);
	if (e_p == 0) return;

	Edge *e = this->OnConvexHull(p);
	if (e == 0)
	{
		e = e_p->dnext;
		this->RemovePoint(p);
		this->TriangulateEmptyPolygon(*e);
	}
	else
	{
		GridPoint end_point = e->twin->dnext->twin->orig;
		Edge *start_edge = e->dnext;
		this->RemovePoint(p);
		this->TriangulateBorder(*start_edge, end_point);
	}
}

template class GridTriangulation<SparseGrid, SparsePointIter>;
template class GridTriangulation<FullGrid, FullPointIter>;
