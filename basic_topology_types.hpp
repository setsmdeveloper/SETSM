#ifndef BASIC_TOPOLOGY_TYPES_H
#define BASIC_TOPOLOGY_TYPES_H

#include <cstddef>
#include <cstdint>
#include <vector>

/**************************
 **************************
 * MACROS, TYPEDEFS,      *
 * FOREWARD DECLARATIONS, *
 * ETC.                   *
 **************************
 **************************/
#define THREADED_CUTOFF 100	// Number of points after which the algorithm should switch to a serial triangulation
#define RETRI_CUTOFF 50	// Number of points after which the algorithm should switch to a serial triangulation
#define TINUPD_THRSHLD 100 // Used to decide update vs. new triangulation; larger means less likely to use TINUpdate

// We don't really need 128-bit integers but here's the code in case it's ever wanted
typedef int64_t INT64;
//#ifdef __SIZEOF_INT128__ // __int128_t is a compiler extension
//	typedef __int128_t INT128; 
//#else // Need to implement own 128-bit int type (or live with 64-bit?)
	typedef int64_t INT128;
//#endif
typedef int INDEX;
typedef struct GridPoint GridPoint;
typedef struct Edge Edge;
typedef struct ExtremeEdges ExtremeEdges;

/*************************************************
 *************************************************
 * Primary data structures used in triangulation *
 *************************************************
 *************************************************/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Data structure for point in grid
struct GridPoint
{
	INDEX row; // Coordinates
	INDEX col;
	bool operator==(const GridPoint &other) { return this->row == other.row && this->col == other.col; };
	bool operator!=(const GridPoint &other) { return this->row != other.row || this->col != other.col; };
};

// Compare function to sort by x/col then y/row
inline bool LessThanXY(const GridPoint &a, const GridPoint &b) { return (a.col < b.col) || ((a.col == b.col) && (a.row < b.row)); }
// Compare function to sort by y/row then x/col
inline bool LessThanYX(const GridPoint &a, const GridPoint &b) { return (a.row < b.row) || ((a.row == b.row) && (a.col > b.col)); } 
// Conversion functions from GridPoint to linear index and vice versa
inline GridPoint Convert(std::size_t index, INDEX width) { GridPoint p = { (int)(index / width), (int)(index % width) }; return p; }
inline std::size_t Convert(const GridPoint &p, INDEX width) { return p.row * width + p.col; }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Data structure for (half-)edge in grid
struct Edge
{
	Edge *dnext;	// Next edge out of destination vertex, counter-clockwise. That is, starting from a point on Edge e, rotate CCW. The next edge you tough this way is e's dnext. Note that this is opposite from some sources where dnext would be the next clockwise edge.
	Edge *oprev;	// Previous edge into origin vertex, counter-clockwise (that is, next edge into origin, clockwise)
	Edge *twin;	// The companion edge with origin and destination switched
	GridPoint orig;	// The origin of the edge
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Data structure to hold various extreme edges of a triangulation/convex hull
struct ExtremeEdges
{
	Edge *left_edge_ccw;	// The counter-clockwise edge of the hull with leftmost origin
	Edge *right_edge_cw;	// The clockwise edge of the hull with rightmost origin
	Edge *bottom_edge_ccw;	// The counter-clockwise edge of the hull with lowest origin
	Edge *top_edge_cw;	// The clockwise edge of the hull with highest origin
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Forward declaration
template <typename GridType, typename IterType>
class Grid;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Data structure to store the (half-)edges of a triangulation
class EdgeList
{
	private:
		Edge *edges;		// List of Edge objects
		// List of pointers to currently unused positions in 'edges'
		union {
			struct { Edge **unused_edges; std::size_t idx; };
			std::vector<Edge *> *local_copy_unused_edges;
		};
		std::size_t size;	// Maximum number of edges to store
		bool is_part_of_split;	// Denotes whether this is part of a larger, split EdgeList
		bool is_local_copy;	// Denotes whether this is a copy of the global EdgeList

		// Alternative constructor to create EdgeList from components of a base EdgeList
		//   edges - list of edge objects from base EdgeList
		//   unused_edges - list of unused edges from base EdgeList
		//   num_points = maximum number of points for THIS triangulation (not for base triangulation)
		EdgeList(Edge *edges, Edge **unused_edges, std::size_t num_points);
		// Alternative constructor to create a local EdgeList from a global EdgeList
		// Used for parallel point removal
		//   global_list - global EdgeList
		//   num_unused = maximum number of unused edges for local point removal
		EdgeList(Edge *edges, std::size_t num_points);
	public:
		// Standard constructor to create EdgeList from scratch
		//   num_points - maximum number of points in corresponding triangulation
		EdgeList(std::size_t num_points);
		// Deconstructor
		~EdgeList();
		// Returns a pointer to a valid Edge object for use in the triangulation
		Edge *GetNewEdge();
		// Removes an Edge from triangulation
		//   edge - the Edge object to be removed from triangulation
		void RemoveEdge(Edge &edge);
		// Splits the current EdgeList into two sub lists
		//   left_list - pointer to an EdgeList; will be replaced with 'left half'
		//     of the current EdgeList
		//   right_list - pointer to an EdgeList; will be replaced with 'right half'
		//     of the current EdgeList
		//   median - the desired number of points in the triangulation corresponding
		//     to left_list (leaving the rest for right_list)
		// The current edge list should be empty (all edges unused) when split for the
		// SplitEdgeList + MergeEdgeLists combination to work properly
		void SplitEdgeList(EdgeList *&left_list, EdgeList *&right_list, std::size_t median);
		// Merges two EdgeLists into one EdgeList (this)
		//   left_list - the left EdgeList from SplitEdgeList call
		//   right_list - the right EdgeList from SplitEdgeList call
		// The input EdgeLists should have been created by calling SplitEdgeList on this
		void MergeEdgeLists(EdgeList &left_list, EdgeList &right_list);
		// Creates a copy of the current EdgeList with its own empty unused_edges and idx
		EdgeList *MakeLocalCopy();
		// Merges local unused list into the current unused list
		//   local_list - the left EdgeList from MakeLocalCopy call
		void MergeUnused(EdgeList &local_list);

		// Sets grid so that each grid point corresponding to a point in the
		// triangulation stores an pointer to an Edge object out of that point
		//   grid - the grid object for storing Edge pointers
		// Best to call only once after EdgeList (corresponding GridTriangulation)
		// has been fully triangulated
		template <typename GridType, typename IterType>
		void SetGrid(Grid<GridType, IterType> *grid);
};

template <typename GridType, typename IterType>
void EdgeList::SetGrid(Grid<GridType, IterType> *grid)
{
	bool *invalid_edge = new bool[this->size];
	#pragma omp parallel
	{
		// Create a companion array that holds
		// booleans tracking whether or not the
		// corresponding edge memory is actually
		// used
		#pragma omp for
		for (size_t t = 0; t < this->size; t++)
		{
			invalid_edge[t] = false;
		}
		#pragma omp for
		for (size_t t = 0; t < this->idx; t++)
		{
			invalid_edge[this->unused_edges[t] - this->edges] = true;
		}

		// Loop over edges and set the edge out of the
		// origin of the current edge to the current edge
		// Threading race conditions are not an issue, since
		// and succesful write is valid (TODO ?)
		/*
		#pragma omp for
		for (size_t t = 0; t < this->size; t++)
		{
			if (invalid_edge[t]) continue;
			
			grid->SetElem((this->edges[t]).orig, this->edges + t);
		}
		*/
	}

	for (size_t t = 0; t < this->size; t++)
	{
		if (invalid_edge[t]) continue;

		
		grid->SetElem((this->edges[t]).orig, this->edges + t);
	}

	delete [] invalid_edge;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // BASIC_TOPOLOGY_TYPES_H
