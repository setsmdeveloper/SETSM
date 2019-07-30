#ifndef GRID_TRIANGULATION_H
#define GRID_TRIANGULATION_H

#include <cstddef>
#include <stdint.h>
#include <vector>

#include "basic_topology_types.hpp"
#include "grid_iterators.hpp"
#include "grid_types.hpp"

template <typename GridType, typename IterType>
class GridTriangulation
{
	private:
		Grid<GridType, IterType> *grid;
		EdgeList *edge_list;
		INDEX width;
		INDEX height;

		Edge *AddEdgeAndTwin(const GridPoint &orig, const GridPoint &dest);
		void RemoveEdgeAndTwin(Edge &edge);
		void Weld(Edge &in, Edge &out);
		Edge *Bridge(Edge &start, Edge &end);	
		ExtremeEdges *TriangulateHorizontalSerial(GridPoint *points[], size_t num_points);
		ExtremeEdges *TriangulateHorizontalThreaded(GridPoint *points[], size_t num_points);
		ExtremeEdges *TriangulateVerticalSerial(GridPoint *points[], size_t num_points);
		ExtremeEdges *TriangulateVerticalThreaded(GridPoint *points[], size_t num_points);
		ExtremeEdges *Triangulate2(GridPoint *points[]);
		ExtremeEdges *Triangulate3(GridPoint *points[]);
		Edge *NextCrossEdge(Edge *base);
		Edge *OnConvexHull(const GridPoint &p);
		
	public:
		GridTriangulation(INDEX width, INDEX height);
		~GridTriangulation();
		Grid<GridType, IterType> *getGrid();
		void Triangulate(GridPoint *points[], size_t num_points);
		Edge *FindFaceRepresentative(Edge *edge);
		void RemoveAllPoints(GridPoint **blunders, size_t num_blunders);
		void Retriangulate(GridPoint **blunders, size_t num_blunders);
		void RetriangulateVertical(std::vector<GridPoint*> blunders, size_t num_blunders);
		void RetriangulateHorizontal(std::vector<GridPoint*> blunders, size_t num_blunders);
		void SetEdgeOut(const GridPoint &p, Edge *e) { this->grid->SetElem(p, e); }
		Edge *GetEdgeOut(const GridPoint &p) { return this->grid->GetElem(p); }
		void RemovePoint(const GridPoint &p)
		{
			Edge *e = this->GetEdgeOut(p);
			if (e == 0) return;

			std::vector<Edge *> edges;
			Edge *f = e;
			do { edges.push_back(f); } while ((f = f->twin->dnext) != e);
			for (auto it = edges.begin(); it != edges.end(); ++it)
			{
				Edge *e = *it;
				this->SetEdgeOut(p, ((e->twin->dnext == e) ? 0 : e->twin->dnext));
				this->SetEdgeOut(e->twin->orig, ((e->dnext == e->twin) ? 0 : e->dnext));
				this->RemoveEdgeAndTwin(**it);
			}
			this->grid->IgnorePoint(p);
		}
		void TriangulateEmptyPolygon(Edge &edge);
		void TriangulateEmpty4(Edge &edge);	
		void TriangulateEmpty5(Edge &edge);	
		void TriangulateEmpty6(Edge &edge);
		void TriangulateEmpty7(Edge &edge);
		void TriangulateBorder(Edge& e, const GridPoint& p);
		void RemovePointAndRetriangulate(const GridPoint &p);

		std::size_t GetAllTris(GridPoint tris[][3])
		{
			if (this->grid->NumOfPts() < 3)
			{
				return 0;
			}
			else
			{
				return this->grid->GetAllTris(tris);
			}
		}
		std::size_t NumOfPts() { return this->grid->NumOfPts(); }
};

typedef GridTriangulation<SparseGrid, SparsePointIter> SparseTriangulation;
typedef GridTriangulation<FullGrid, FullPointIter> FullTriangulation;

template void EdgeList::SetGrid(Grid<SparseGrid, SparsePointIter> *grid);
template void EdgeList::SetGrid(Grid<FullGrid, FullPointIter> *grid);

#endif	// GRID_TRIANGULATION_H
