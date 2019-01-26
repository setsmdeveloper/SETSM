#ifndef GRID_TYPES_H
#define GRID_TYPES_H

#include "basic_topology_types.hpp"
#include "grid_iterators.hpp"

// Grid type class to represent 'full' grid (full matrix representation)
class FullGrid
{
	private:
		Edge **elems;		// Grid storage (list of pointers to Edges)
		std::size_t num_points;	// Current number of set points in grid
		const INDEX width;	// Grid dimensions
		const INDEX height;
	public:
		// Standard constructor
		//   width - width of grid
		//   height - height of grid
		//   points - list of pointers to points in grid to be triangulated (only here for consistency with SparseGrid constructor)
		//   num_points - number of points in grid to be triangulated
		FullGrid(INDEX width, INDEX height, GridPoint *points[], std::size_t num_points):width(width),height(height),num_points(num_points) { this->elems = new Edge*[width * height](); }
		~FullGrid() { delete [] this->elems; }
		// Sets the the edge out of point in grid
		//   p - point in grid
		//   e - edge out of p
		void SetElem(const GridPoint& p, Edge *e) { this->elems[Convert(p, this->width)] = e; }
		// Retrieves the edge out of point in grid
		//   p - point to retrieve edge out of
		Edge *GetElem(const GridPoint& p) { return this->elems[Convert(p, this->width)]; }
		// No longer consider point in triangulation
		//   p - point to be removed from grid's consideration
		void IgnorePoint(const GridPoint& p) { this->num_points--; this->SetElem(p, 0); }
		std::size_t NumOfElems() { return this->num_points; }

		// Start and finish iterators over the points in the grid
		FullPointIter PointBegin() { return FullPointIter(this->elems, this->width, this->height, 0); }
		FullPointIter PointEnd() { return FullPointIter(this->elems, this->width, this->height, this->width * this->height); }
};

// Grid type class to represent 'sparse' grid (hashtable representation)
class SparseGrid
{
	private:
		std::unordered_map<std::size_t, std::size_t> coord_2_index;	// Hashtable from index in grid (row*width+col) to index in elems
		Edge **elems;							// Grid storage (list of pointers to Edges)
		INDEX width;							// Grid width
	public:
		// Standard constructor
		//   width - width of grid
		//   height - height of grid (only here for consistency with FullGrid constructor)
		//   points - list of pointers to points in grid to be triangulated
		//   num_points - number of points in grid to be triangulated
		// Underlying hashtable is initialized here and should be unchanged from here on out
		SparseGrid(INDEX width, INDEX height, GridPoint *points[], std::size_t num_points):width(width)
		{
			this->elems = new Edge*[num_points];
			this->coord_2_index.reserve(num_points);
			for (std::size_t t = 0; t < num_points; t++) coord_2_index[Convert(*(points[t]), this->width)] = t;
		}
		~SparseGrid() { delete [] this->elems; }
		// Sets the the edge out of point in grid
		//   p - point in grid
		//   e - edge out of p
		// Only call on points used in construction of this
		void SetElem(const GridPoint &p, Edge *e) { this->elems[coord_2_index.at(Convert(p, this->width))] = e; }
		// Retrieves the edge out of point in grid
		//   p - point to retrieve edge out of
		// Only call on points used in construction of this
		Edge *GetElem(const GridPoint &p) { return this->elems[coord_2_index.at(Convert(p, this->width))]; }
		// No longer consider point in triangulation
		//   p - point to be removed from grid's consideration
		void IgnorePoint(const GridPoint &p) { this->coord_2_index.erase(Convert(p, this->width)); }
		std::size_t NumOfElems() { return this->coord_2_index.size(); }

		// Start and finish iterators over the points in the grid
		SparsePointIter PointBegin() { return SparsePointIter(this->width, coord_2_index.begin()); }
		SparsePointIter PointEnd() { return SparsePointIter(this->width, coord_2_index.end()); }
};

// Grid class for storing edges out of points in triangulation
template <typename GridType, typename IterType>
class Grid
{
	private:
		GridType *grid;	// Underlying storage container
	public:
		// Standard constructor
		//   width - width of grid
		//   height - height of grid
		//   points - list of pointers to points in grid to be triangulated
		//   num_points - number of points in grid to be triangulated
		Grid(INDEX width, INDEX height, GridPoint *points[], std::size_t num_points);
		~Grid();
		// Sets the the edge out of point in grid
		//   p - point in grid
		//   e - edge out of p
		// Only call on points used in construction of this
		void SetElem(const GridPoint &p, Edge *e) { this->grid->SetElem(p, e); }
		// Retrieves the edge out of point in grid
		//   p - point to retrieve edge out of
		// Only call on points used in construction of this
		Edge *GetElem(const GridPoint &p) { return this->grid->GetElem(p); }
		// No longer consider point in triangulation
		//   p - point to be removed from grid's consideration
		void IgnorePoint(const GridPoint &p) { this->grid->IgnorePoint(p); }
		// Fills list with triangles of current triangulation
		//   tris - already created list to store tris. should be
		//     at least 2*num_points long
		std::size_t GetAllTris(GridPoint tris[][3]);
		std::size_t NumOfPts() { return this->grid->NumOfElems(); }

		// Iterators over both
		//   a) points currently in triangulation
		//   b) triangles adjacent to a given point in current triangulation
		IterType PointBegin() { return this->grid->PointBegin(); }
		IterType PointEnd() { return this->grid->PointEnd(); }
		AdjTriIter AdjTriBegin(const GridPoint &p) { Edge *e = this->GetElem(p); return AdjTriIter(true, e, e); }
		AdjTriIter AdjTriEnd(const GridPoint &p) { Edge *e = this->GetElem(p); return AdjTriIter(false, e, e); }
};

#endif // GRID_TYPES_H
