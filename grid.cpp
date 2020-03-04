#include "grid_types.hpp"

template <typename GridType, typename IterType>
Grid<GridType, IterType>::Grid(INDEX width, INDEX height, GridPoint *points[], std::size_t num_points)
{
	this->grid = new GridType(width, height, points, num_points);
}

template <typename GridType, typename IterType>
Grid<GridType, IterType>::~Grid()
{
	delete this->grid;
}

template <typename GridType, typename IterType>
std::size_t Grid<GridType, IterType>::GetAllTris(std::vector<Tri> *tris)
{
	// Loop over all points in triangulation, collecting adjacent triangles
	// To prevent 'over-counting' only collect triangle from its smallest vertex
	// Should be able to be parallelized, but haven't yet found a good one. 
	// Be wary, list order does matter
	for (auto p_it = this->PointBegin(); p_it != this->PointEnd(); p_it++)
	{
		Edge *e = this->GetElem(*p_it);
		Edge *f = e->twin->dnext->twin->dnext; // If degree of current point is 1 or 2, then f == e

		GridPoint p1 = e->orig;
		GridPoint p2 = e->twin->orig;
		GridPoint p3 = e->twin->dnext->twin->orig;

		// Check that degree > 2 or
		//  at least edge(s) out not parallel
		if ((f != e) || !((p1.col - p2.col)*(p1.row - p3.row) == (p1.row - p2.row)*(p1.col - p3.col)))
		{
			for (auto t_it = this->AdjTriBegin(*p_it); t_it != this->AdjTriEnd(*p_it); t_it++)
			{
				GridPoint *points = *t_it;
				if (LessThanXY(points[0], points[1]) && LessThanXY(points[0], points[2]))
				{
					Tri temp_tri;
					temp_tri.pts[0] = points[0];
					temp_tri.pts[1] = points[1];
					temp_tri.pts[2] = points[2];
					tris->push_back(temp_tri);
				}
			}
		}
	}

	return tris->size();
}

template class Grid<SparseGrid, SparsePointIter>;
template class Grid<FullGrid, FullPointIter>;
