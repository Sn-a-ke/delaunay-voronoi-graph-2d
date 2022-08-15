// Copyright 2020 Alexandr Marchenko. All Rights Reserved.

#pragma once

#include <vector>
#include <algorithm>
#include "vector_graph.h"
#include "vector_graph_delaunay.h"


namespace helper
{
	template<typename ElementType>
	static uint32_t find_index(const std::vector<ElementType>& Arr, const ElementType& Element)
	{
		auto f = std::find(Arr.begin(), Arr.end(), Element);
		if (f != Arr.end())
		{
			return static_cast<uint32_t>(f - Arr.begin());
		}
		return -1;
	}

	template<typename ElementType>
	static uint32_t add_unique(std::vector<ElementType>& Arr, ElementType&& Element)
	{
		const uint32_t Idx = find_index(Arr, Element);
		if (Idx != -1)
		{
			return Idx;
		}
		Arr.push_back(std::move(Element));
		return Arr.size() - 1;
	}

	template<typename Real, typename VectorType>
	static VectorType find_circum_center_2d(const VectorType& A, const VectorType& B, const VectorType& C)
	{
		const Real Determinant = 2.0 * (A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]));
		if (Determinant != static_cast<Real>(0.0))
		{
			const Real Asq = A[0] * A[0] + A[1] * A[1];//size_squared_2d(A);
			const Real Bsq = B[0] * B[0] + B[1] * B[1];//size_squared_2d(B);
			const Real Csq = C[0] * C[0] + C[1] * C[1];//size_squared_2d(C);

			const Real Ux = (Asq * (B[1] - C[1]) + Bsq * (C[1] - A[1]) + Csq * (A[1] - B[1])) / Determinant;
			const Real Uy = (Asq * (C[0] - B[0]) + Bsq * (A[0] - C[0]) + Csq * (B[0] - A[0])) / Determinant;

			VectorType Out = VectorType();
			Out[0] = Ux;
			Out[1] = Uy;
			return Out;
		}
		__assume(0);
		return VectorType();
	}

	template<typename VectorType>
	static VectorType A_to_B_left_ray(const VectorType& A, const VectorType& B, double Length = 1.0)
	{
		const VectorType V = B - A;
		VectorType Vd = VectorType();
		Vd[0] = V[1];
		Vd[1] = -V[0];
		const VectorType Mid = (V * 0.5) + A;
		const double Scale = Length / std::sqrt(Vd[0] * Vd[0] + Vd[1] * Vd[1]); // todo:  fast InvSqrt
		return Mid + (Vd * Scale);
	}

} // namespace helper

/** Voronoi */ //todo: voronoi_graph constructor from points
template<typename InVectorType, typename IntType = uint16_t>
class voronoi_graph
{
public:
	using VectorType = InVectorType; //todo

	using VectorGraphType = TVectorOrientedGraph<VectorType, IntType, std::vector<IntType>>;
	using VoronoiGraphType = TVectorOrientedGraph<VectorType, IntType, std::vector<IntType>>;

	using VoronoiLinkType = std::vector<IntType>;
	using VoronoiLinkContainer = std::vector<VoronoiLinkType>;


	explicit voronoi_graph<InVectorType, IntType>(VectorGraphType& InOriginDelaunay, const bool bInit = true)
		: OriginDelaunay(InOriginDelaunay)
		, VoronoiGraph(VoronoiVerts)
	{
		if (bInit)
		{
			BuildVoronoi_Unsorted();
			//BuildVoronoi_V2();
		}
	}

	/*****************************/

	VectorGraphType& OriginDelaunay;

	std::vector<VectorType> VoronoiVerts;
	VoronoiGraphType VoronoiGraph;
	std::vector<std::vector<IntType>> OriginToVoronoi; // todo:: OriginToVoronoi inline allocator

	/*****************************/

	//sorted
	void BuildVoronoi_V2(/*FBox2D* InClipBox = nullptr*/) //todo: custom clip shape
	{
		const int32_t OriginNum = OriginDelaunay.Nodes.size();

		VoronoiGraph.NodeLinks.reserve(OriginNum * 2 + 10);
		VoronoiVerts.reserve(OriginNum * 2 + 2);

		OriginToVoronoi.resize(OriginNum, std::vector<IntType>());
		for (int32_t i = 0; i < OriginNum; ++i)
		{
			OriginToVoronoi[i].resize(OriginDelaunay.NodeLinks[i].size(), std::numeric_limits<IntType>::max());
		}

		for (int32_t i = 0; i < OriginNum; ++i)
		{
			const VectorType& V0 = OriginDelaunay.Nodes[i];

			const auto& Lin = OriginDelaunay.NodeLinks[i];
			const int32_t Num = Lin.size();

			// iterate links for current_link_Id and next_link_Id
			const bool bBorder = OriginDelaunay.IsBorder(i);
			int32_t N2 = static_cast<int32_t>(bBorder);
			int32_t N1 = helper::prev_index(N2, Num);
			for (; N2 < Num; N1 = N2++)
			{
				const IntType L1 = Lin[N1];						 // link 1
				const IntType L2 = Lin[N2];						 // link 2
				const VectorType& V1 = OriginDelaunay.Nodes[L1]; // vert 1
				const VectorType& V2 = OriginDelaunay.Nodes[L2]; // vert 2
				if (L1 > i && L2 > i)							 // link 1 2 > i
				{
					VoronoiVerts.push_back(helper::find_circum_center_2d(V0, V1, V2)); // add voronoi vert
					VoronoiGraph.NodeLinks.push_back(std::vector<IntType>());		   // add voronoi link array
					const IntType VoroNodeIdx = VoronoiVerts.size() - 1;

					const auto& Lk_1 = OriginDelaunay.NodeLinks[L1]; // link container for link 1
					const auto& Lk_2 = OriginDelaunay.NodeLinks[L2]; // link container for link 2

					const int32_t Lid1 = helper::find_index(Lk_1, L2);
					const int32_t Lid2 = helper::prev_index(helper::find_index(Lk_2, L1), Lk_2.Num());

					OriginToVoronoi[i][N1] = VoroNodeIdx;
					OriginToVoronoi[L1][Lid1] = VoroNodeIdx;
					OriginToVoronoi[L2][Lid2] = VoroNodeIdx;
				}
			}
		}

		for (int32_t i = 0; i < OriginNum; ++i)
		{
			const auto& OrToVoro_i = OriginToVoronoi[i];
			const bool bBorder = OrToVoro_i.back() == std::numeric_limits<IntType>::max();
			const int32_t Num2 = OrToVoro_i.size() - static_cast<int32_t>(bBorder);
			for (int32_t j = 0; j < Num2; ++j)
			{
				const IntType V_Idx_1 = OrToVoro_i[j];
				const IntType V_Idx_2 = OrToVoro_i[helper::next_index(j, Num2)];
				helper::add_unique(VoronoiGraph.NodeLinks[V_Idx_1], V_Idx_2); //todo BuildVoronoi VoronoiGraph.NodeLinks AddUnique optimize
				helper::add_unique(VoronoiGraph.NodeLinks[V_Idx_2], V_Idx_1);
			}
			//VoronoiGraph.SortAdjacent(i);
		}
	}

	void BuildVoronoi_Unsorted(const double Offset = 1000.0 /*FBox2D* InClipBox = nullptr*/) //todo: custom clip shape
	{
		const int32_t OriginNum = OriginDelaunay.Nodes.size();

		std::vector<IntType> SortId; // = OriginDelaunay.SortedArrayView2d();
		SortId.reserve(OriginDelaunay.Nodes.size());

		OriginToVoronoi.resize(OriginNum, std::vector<IntType>());

		for (int32_t i = 0; i < OriginNum; ++i)
		{
			SortId.push_back(i);
			const bool bBorder = OriginDelaunay.IsBorder(Idx);
			const int32_t Size = OriginDelaunay.NodeLinks[i].size() + static_cast<int32_t>(bBorder);
			OriginToVoronoi[i].resize(Size, std::numeric_limits<IntType>::max());
		}

		std::sort(
			SortId,
			[&](const IntType& A, const IntType& B)
			{
				const VectorType& V1 = OriginDelaunay.Nodes[A];
				const VectorType& V2 = OriginDelaunay.Nodes[B];
				return V1[0] < V2[0] || (V1[0] == V2[0] && V1[1] < V2[1]);
			});

		VoronoiVerts.reserve(OriginNum * 2 + 2);
		VoronoiGraph.NodeLinks.reserve(OriginNum * 2 + 10);

		for (const int32_t Idx : SortId)
		{
			const VectorType& V0 = OriginDelaunay.Nodes[Idx];
			const auto& CurrentLinks = OriginDelaunay.NodeLinks[Idx];

			const bool bBorder = CurrentLinks.size() != OriginToVoronoi[Idx].size();
			if (bBorder)
			{
				auto& OrToVorIdx = OriginToVoronoi[Idx];

				const IntType L1 = CurrentLinks[0];
				if (L1 > Idx)
				{
					const VectorType& V1 = OriginDelaunay.Nodes[L1];
					VectorType P1 = helper::A_to_B_left_ray(V0, V1, Offset);

					VoronoiVerts.push_back(std::move(P1));
					VoronoiGraph.NodeLinks.push_back(std::vector<IntType>());
					const IntType VoroNodeIdx = = VoronoiVerts.size() - 1;

					auto& OrToVorL1 = OriginToVoronoi[L1];

					OrToVorIdx.back() = VoroNodeIdx;
					OrToVorL1[OrToVorL1.size() - 2] = VoroNodeIdx;
				}

				const IntType L2 = CurrentLinks.back();
				if (L2 > Idx)
				{
					const VectorType& V2 = OriginDelaunay.Nodes[L2];
					VectorType P2 = helper::A_to_B_left_ray(V2, V0, Offset);

					VoronoiVerts.push_back(std::move(P2));
					VoronoiGraph.NodeLinks.push_back(std::vector<IntType>());
					const IntType VoroNodeIdx = = VoronoiVerts.size() - 1;

					auto& OrToVorL2 = OriginToVoronoi[L2];

					OrToVorL2.back() = VoroNodeIdx;
					OrToVorIdx[OrToVorIdx.size() - 2] = VoroNodeIdx;
				}
			}

			// circum center for each triangle
			const int32_t Size = CurrentLinks.size();
			int32_t N2 = static_cast<int32_t>(bBorder);
			int32_t N1 = helper::prev_index(N2, Size);
			for (; N2 < Size; N1 = N2++)
			{
				const IntType L1 = CurrentLinks[N1];
				const IntType L2 = CurrentLinks[N2];
				if (L1 > Idx && L2 > Idx)
				{
					const VectorType& V1 = OriginDelaunay.Nodes[L1];
					const VectorType& V2 = OriginDelaunay.Nodes[L2];

					VoronoiVerts.push_back(helper::find_circum_center_2d(V0, V1, V2));
					VoronoiGraph.NodeLinks.push_back(std::vector<IntType>());

					const IntType VoroNodeIdx = VoronoiVerts.size() - 1;

					const auto& Lk_1 = OriginDelaunay.NodeLinks[L1];
					const auto& Lk_2 = OriginDelaunay.NodeLinks[L2];

					const int32_t Lid1 = helper::find_index(Lk_1, L2);
					const int32_t Lid2 = helper::prev_index(helper::find_index(Lk_2, L1), Lk_2.size());

					OriginToVoronoi[Idx][N1] = VoroNodeIdx;
					OriginToVoronoi[L1][Lid1] = VoroNodeIdx;
					OriginToVoronoi[L2][Lid2] = VoroNodeIdx;
				}
			}
		}

		for (const int32_t Idx : SortId)
		{
			const auto& OrToVoro_Idx = OriginToVoronoi[Idx];
			const int32_t Size = OrToVoro_Idx.size();
			for (int32_t j = 0; j < Size; j++)
			{
				const IntType V_Idx_1 = OrToVoro_Idx[j];
				const IntType V_Idx_2 = OrToVoro_Idx[helper::next_index(j, Size)];
				helper::add_unique(VoronoiGraph.NodeLinks[V_Idx_1], V_Idx_2); //todo BuildVoronoi VoronoiGraph.NodeLinks AddUnique optimize
				helper::add_unique(VoronoiGraph.NodeLinks[V_Idx_2], V_Idx_1);
			}
		}
	}

	/**
	* Lloyd's algorithm Voronoi iteration or relaxation
	* https://en.wikipedia.org/wiki/Lloyd%27s_algorithm
	*/
	std::vector<VectorType> VoronoiRelax()
	{
		/*
		VectorType Centroid() const
		{
			constexpr float f = 1.0 / 3.0;
			return VectorType(
				(V[0].X + V[1].X + V[2].X) * f,
				(V[0].Y + V[1].Y + V[2].Y) * f,
				(V[0].Z + V[1].Z + V[2].Z) * f
			);
		}
		*/

		const int32_t Num = OriginDelaunay.Nodes.size();

		std::vector<VectorType> Out; //todo self relax?
		Out.resize(Num, VectorType());

		for (int32_t i = 0; i < Num; ++i)
		{
			float Area = 0;
			const VectorType& A = OriginDelaunay.Nodes[i];
			auto& VoroIds = OriginToVoronoi[i];
			std::vector<std::tuple<float, VectorType>> SSCen;
			SSCen.reserve(VoroIds.size());

			for (int32_t j = 0; j < VoroIds.size(); ++j)
			{
				const int32_t Idx_1 = VoroIds[j];
				const int32_t Idx_2 = VoroIds[GeomHelpers::next_index(j, VoroIds.size())];
				const VectorType& B = VoronoiVerts[Idx_1];
				const VectorType& C = VoronoiVerts[Idx_2];

				const float TriArea = std::abs(A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1])) * 0.5f;
				Area += TriArea;

				SSCen.push_back(std::make_tuple(TriArea, (A + B + C) / 3.f));
			}

			for (int32_t j = 0; j < VoroIds.size(); ++j)
			{
				const std::tuple<float, VectorType>& Tup = SSCen[j];
				Out[i] += (std::get<1>(Tup) * std::get<0>(Tup)) / Area;
			}
		}
		return Out;
	}
};


namespace Triangulate
{
	template<typename VectorType, typename IntType = uint16_t, typename LinkType = std::vector<IntType>>
	static voronoi_graph<VectorType, IntType> triangulate(TVectorOrientedGraph<VectorType, IntType, LinkType>& TriangulatedGraph)
	{
		return voronoi_graph<VectorType, IntType>(TriangulatedGraph, true);
	}
} // namespace Triangulate