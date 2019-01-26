#include <cstddef>
#include "basic_topology_types.hpp"

EdgeList::EdgeList(std::size_t num_points)
{
	// 3F < 2E, E <= F + V <= (2/3)E + V <-> E <= 3V
	// We store edge and twin, so need memory of 6*num_points
	// Allocate all memory upfront and use stack structure
	// (this->unused_edges + this->idx) to track currently
	// unused memory
	this->size = 6 * num_points;
	this->idx = this->size;
	this->edges = new Edge[this->size];
	this->unused_edges = new Edge*[this->size];
	for (std::size_t t = 0; t < this->size; t++)
	{
		this->unused_edges[t] = this->edges + t;
	}
	this->is_part_of_split = false;
}

EdgeList::EdgeList(Edge *edges, Edge **unused_edges, size_t num_points)
{
	// Constructor used to make sub-EdgeList
	// from a larger EdgeList
	this->size = 6 * num_points;
	this->idx = this->size;
	this->edges = edges;
	this->unused_edges = unused_edges;
	this->is_part_of_split = true;
}

EdgeList::~EdgeList()
{
	// In case this is not part of a
	// split EdgeList, free edge memory
	if (!this->is_part_of_split)
	{
		delete [] this->edges;
		delete [] this->unused_edges;
	}
}

Edge *EdgeList::GetNewEdge()
{
	this->idx--;
	return this->unused_edges[this->idx];
}

void EdgeList::RemoveEdge(Edge &edge)
{
	this->unused_edges[this->idx] = &edge;
	this->idx++;
}

void EdgeList::SplitEdgeList(EdgeList *&left_list, EdgeList *&right_list, size_t median)
{
	// Splits this into two sub EdgeLists with enough memory to store
	// edges for a triangulation of median and (this->size / 6) - median
	// points, respectively
	left_list = new EdgeList(this->edges, this->unused_edges, median);
	right_list = new EdgeList(this->edges + 6 * median, this->unused_edges + 6 * median, (this->size / 6) - median);
}

void EdgeList::MergeEdgeLists(EdgeList &left_list, EdgeList &right_list)
{
	// Combine lists of unused edges, updating
	// this->idx appropriately
	this->idx = left_list.idx;
	for (size_t t = 0; t < right_list.idx; t++)
	{
		this->unused_edges[this->idx + t] = right_list.unused_edges[t];
	}
	this->idx += right_list.idx;
}
