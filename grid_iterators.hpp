#ifndef GRID_ITERATORS_H
#define GRID_ITERATORS_H

#include "basic_topology_types.hpp"
#include <unordered_map>

// Iterator class to iterate over full grid (full matrix representation)
class FullPointIter : public std::iterator<std::forward_iterator_tag, GridPoint>
{
	private:
		Edge **elems;		// List of Edge pointers
		const INDEX width;	// Grid dimensions
		const INDEX height;
		const std::size_t max;	// Length of grid = width * height for convenience
		std::size_t curr;	// Current index in grid

		// Increment grid position while current element is null until
		//   a) hit end of grid
		//   b) current element is no longer null
		void UpdateUntilValid() { while ((this->curr < this->max) && (this->elems[this->curr] == 0)) this->curr++; }
	public:
		// Standard constructor
		//   elems - underlying list of Edge pointers from grid
		//   width - width of grid
		//   height - height of grid
		//   curr - starting index in grid (0 = begin, width*height = end) 
		FullPointIter(Edge **elems, INDEX width, INDEX height, std::size_t curr):elems(elems),width(width),height(height),curr(curr),max(width*height) { this->UpdateUntilValid(); }
		~FullPointIter() {}
		bool operator==(const FullPointIter &other) { return (curr == other.curr) && (width == other.width) && (height == other.height) && (elems == other.elems); }
		bool operator!=(const FullPointIter &other) { return (curr != other.curr) || (width != other.width) || (height != other.height) || (elems != other.elems); }
		GridPoint operator*() { return Convert(this->curr, this->width); }
		FullPointIter &operator++() { this->curr++; this->UpdateUntilValid(); return *this; }
		FullPointIter operator++(int) { FullPointIter temp(*this); this->curr++; this->UpdateUntilValid(); return temp; }
};

// Iterator class to iterate over sparse grid (sparse matrix representation)
class SparsePointIter : public std::iterator<std::forward_iterator_tag, GridPoint>
{
	private:
		std::unordered_map<std::size_t, std::size_t>::iterator it;	// Iterator over underlying hashtable of sparse grid
		const INDEX width;						// Width of grid
	public:
		// Stardard constructor
		//   it - iterator over underlying hashtable of grid
		//   width - width of grid
		SparsePointIter(INDEX width, std::unordered_map<std::size_t, std::size_t>::iterator it):width(width),it(it) {}
		~SparsePointIter() {}
		bool operator==(const SparsePointIter &other) { return (width == other.width) && (it == other.it); }
		bool operator!=(const SparsePointIter &other) { return (width != other.width) || (it != other.it); }
		GridPoint operator*() { return Convert(this->it->first, this->width); }
		SparsePointIter &operator++() { ++this->it; return *this; }
		SparsePointIter operator++(int) { SparsePointIter temp(*this); ++this->it; return temp; }
};

// Iterator class to iterate over triangles adjacent to a particular point in a triangulation
class AdjTriIter : public std::iterator<std::forward_iterator_tag, GridPoint*>
{
	private:
		const Edge *e_init;	// Initial edge out of the point
		Edge *e_curr;		// Current edge out of the point
		bool is_first_edge;	// Track whether still on the first edge or not
		GridPoint points[3];	// List to store triple corresponding to triangle TODO should be vector

		// Advance current edge (counterclockwise) to next edge in ring of edges out of point
		void AdvanceCurrent() { this->e_curr = this->e_curr->twin->dnext; }
		// Addvance current edge until it is first in a triangle (prevent degenerate case of edges on convex hull)
		void UpdateUntilValid() { while (this->e_curr->dnext->dnext->dnext != this->e_curr) this->AdvanceCurrent(); }
	public:
		AdjTriIter(bool is_first, Edge *e_init, Edge *e_curr):is_first_edge(is_first),e_init(e_init),e_curr(e_curr) { this->UpdateUntilValid(); }
		~AdjTriIter() {}
		bool operator==(const AdjTriIter &other) { return (is_first_edge == other.is_first_edge) && (e_init == other.e_init) && (e_curr == other.e_curr); }
		bool operator!=(const AdjTriIter &other) { return (is_first_edge != other.is_first_edge) || (e_init != other.e_init) || (e_curr != other.e_curr); }
		GridPoint *operator*()
		{
			this->points[0] = this->e_curr->orig;
			this->points[1] = this->e_curr->dnext->orig;
			this->points[2] = this->e_curr->dnext->dnext->orig;
			return this->points;
		}
		AdjTriIter &operator++()
		{
			this->is_first_edge = false;
			this->AdvanceCurrent();
			this->UpdateUntilValid();
			return *this;
		}
		AdjTriIter operator++(int)
		{
			AdjTriIter temp(*this);
			this->is_first_edge = false;
			this->AdvanceCurrent();
			this->UpdateUntilValid();
			return temp;
		}	
};

#endif // GRID_ITERATORS_H
