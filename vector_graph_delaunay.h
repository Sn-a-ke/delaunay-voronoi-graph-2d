// Copyright 2020 Alexandr Marchenko. All Rights Reserved.

#pragma once

#include "vector_graph.h"
#include <vector>
#include <algorithm>
#include <ppl.h>
#define _USE_MATH_DEFINES
#include <cmath>


namespace helper
{
	static int32_t next_index(const int32_t Id, const int32_t Size, const int32_t Plus = 1)
	{
		const int32_t Result = Id + Plus;
		return Result < Size ? Result : Result - Size;
	}
	static int32_t prev_index(const int32_t Id, const int32_t Size, const int32_t Minus = 1)
	{
		const int32_t Result = Id - Minus;
		return Result < 0 ? Size + Result : Result;
	}

	template<typename Real, typename VectorType>
	struct circle2d_rsq
	{
		Real CenterX, CenterY, RadiusSquared;

		circle2d_rsq() //
			: CenterX(0.0)
			, CenterY(0.0)
			, RadiusSquared(0.0)
		{}
		circle2d_rsq(const VectorType& Center, Real Radius) //
			: CenterX(Center[0])
			, CenterY(Center[1])
			, RadiusSquared(Radius * Radius)
		{}
		circle2d_rsq(const VectorType& A, const VectorType& B, const VectorType& C)
		{
			//const Real A[2] = {static_cast<Real>(InA[0]), static_cast<Real>(InA[1])};
			//const Real B[2] = {static_cast<Real>(InB[0]), static_cast<Real>(InB[1])};
			//const Real C[2] = {static_cast<Real>(InC[0]), static_cast<Real>(InC[1])};

			const Real Asq = (A[0] * A[0]) + (A[1] * A[1]);
			const Real Bsq = (B[0] * B[0]) + (B[1] * B[1]);
			const Real Csq = (C[0] * C[0]) + (C[1] * C[1]);

			const Real R4 = C[0] - B[0];
			const Real R5 = A[0] - C[0];
			const Real R6 = B[0] - A[0];

			const Real R1 = B[1] - C[1];
			const Real R2 = C[1] - A[1];
			const Real R3 = A[1] - B[1];

			const Real d = static_cast<Real>(2.0) * (A[0] * R1 + B[0] * R2 + C[0] * R3);
			CenterX = (Asq * R1 + Bsq * R2 + Csq * R3) / d;
			CenterY = (Asq * R4 + Bsq * R5 + Csq * R6) / d;

			const Real X = A[0] - CenterX;
			const Real Y = A[1] - CenterY;
			RadiusSquared = (X * X) + (Y * Y);
		}

		bool IsInCircle(const VectorType& V) const
		{
			//assert(RadiusSquared != 0.0);
			const Real X = static_cast<Real>(V[0]) - CenterX;
			const Real Y = static_cast<Real>(V[1]) - CenterY;
			const Real DS = (X * X) + (Y * Y);
			return DS < RadiusSquared;
		}

		//FString ToString() const { return FString::Printf(TEXT("CenterX=%lf, CenterY=%lf, RadiusSquared =%lf "), CenterX, CenterY, RadiusSquared); }
	};

	template<typename T>
	static int32_t rotate_to_index(T& InVector, int32_t Index_to_ZeroIndex)
	{
		if (Index_to_ZeroIndex > InVector.size() / 2)
		{
			std::rotate(InVector.rbegin(), InVector.rbegin() + InVector.size() - Index_to_ZeroIndex, InVector.rend());
			return InVector.size() - Index_to_ZeroIndex;
		}
		else
		{
			std::rotate(InVector.begin(), InVector.begin() + Index_to_ZeroIndex, InVector.end());
			return -Index_to_ZeroIndex;
		}
	}

} // namespace helper

namespace Triangulate
{
	template<typename VectorType, typename IntType = uint16_t>
	class convex_delaunay_2d
	{
	public:
		convex_delaunay_2d<VectorType, IntType>() : l_border_id(-1), r_border_id(-1), is_line(false) {}
		~convex_delaunay_2d<VectorType, IntType>() = default;

		std::vector<IntType> convex_border_id;
		int32_t l_border_id;
		int32_t r_border_id;
		bool is_line;

		template<typename VectorGraph>
		void convex_union(
			convex_delaunay_2d<VectorType, IntType>&& r_convex,
			const VectorGraph& graph,
			std::vector<IntType>& l_queue,
			std::vector<IntType>& r_queue)
		{
			convex_delaunay_2d<VectorType, IntType>& l_convex = *this;

			if (l_convex.is_line && r_convex.is_line)
			{
				const VectorType& V1 = graph.Nodes[l_convex.convex_border_id[l_convex.l_border_id]];
				const VectorType& V2 = graph.Nodes[l_convex.convex_border_id[l_convex.r_border_id]];
				const VectorType& V3 = graph.Nodes[r_convex.convex_border_id[r_convex.l_border_id]];
				const VectorType& V4 = graph.Nodes[r_convex.convex_border_id[r_convex.r_border_id]];

				if (helper::cross2d<double, VectorType>(V1, V3, V2) == 0.0 //
					&& helper::cross2d<double, VectorType>(V1, V4, V2) == 0.0)
				{
					l_queue.push_back(l_convex.convex_border_id.back());
					r_queue.push_back(r_convex.convex_border_id[0]);

					vector_append(l_convex.convex_border_id, std::move(r_convex.convex_border_id));
					l_convex.r_border_id = l_convex.convex_border_id.size() - 1;
					return;
				}
			}

			const bool is_four_points = l_convex.convex_border_id.size() == 2 && r_convex.convex_border_id.size() == 2; //todo remove is_four_points ?
			const int32_t left_close_to_right = is_four_points ? -1 : l_convex.r_border_id;
			const int32_t right_close_to_left = is_four_points ? -1 : r_convex.l_border_id;

			auto Tangent = convex_delaunay_2d<VectorType, IntType>::find_tangents_border_Id_to_convex( //
				graph,
				l_convex,
				r_convex,
				left_close_to_right,
				right_close_to_left);

			l_queue.reserve(l_convex.convex_border_id.size());
			if (!l_convex.is_line)
			{
				const int32_t L_Num = l_convex.convex_border_id.size();
				const uint16_t L = l_convex.convex_border_id[Tangent.L_Clock_Id];
				int32_t L_Turn = helper::rotate_to_index(l_convex.convex_border_id, Tangent.L_Clock_Id);

				if (L_Turn < 0)
				{
					L_Turn += L_Num;
				}

				Tangent.L_Clock_Id = 0;
				Tangent.L_InvClock_Id = helper::next_index(Tangent.L_InvClock_Id, L_Num, L_Turn);
				l_convex.l_border_id = helper::next_index(l_convex.l_border_id, L_Num, L_Turn);

				int32_t Start = Tangent.L_InvClock_Id;
				if (Tangent.L_Clock_Id == Tangent.L_InvClock_Id)
				{
					l_queue.push_back(l_convex.convex_border_id[Start]);
					Start = helper::next_index(Start, L_Num);
				}
				for (int32_t i = Start; i != Tangent.L_Clock_Id; i = helper::next_index(i, L_Num))
				{
					l_queue.push_back(l_convex.convex_border_id[i]);
				}
				l_queue.push_back(L);

				auto First = l_convex.convex_border_id.begin() + Tangent.L_InvClock_Id + 1;
				auto Last = First + l_convex.convex_border_id.size() - (Tangent.L_InvClock_Id + 1);
				l_convex.convex_border_id.erase(First, Last);
			}
			else
			{
				l_convex.left_is_line(Tangent.L_Clock_Id, Tangent.L_InvClock_Id, l_queue);
			}
			l_queue.shrink_to_fit();

			r_queue.reserve(r_convex.convex_border_id.size());
			if (!r_convex.is_line)
			{
				const int32_t R_Num = r_convex.convex_border_id.size();
				int32_t R_Turn = helper::rotate_to_index(r_convex.convex_border_id, Tangent.R_Clock_Id);

				if (R_Turn < 0)
				{
					R_Turn += R_Num;
				}

				Tangent.R_Clock_Id = 0;
				Tangent.R_InvClock_Id = helper::next_index(Tangent.R_InvClock_Id, R_Num, R_Turn);
				r_convex.r_border_id = helper::next_index(r_convex.r_border_id, R_Num, R_Turn);

				int32_t Start = Tangent.R_Clock_Id;
				if (Tangent.R_Clock_Id == Tangent.R_InvClock_Id)
				{
					r_queue.push_back(r_convex.convex_border_id[Start]);
					Start = helper::prev_index(Start, R_Num);
				}
				for (int32_t i = Start; i != Tangent.R_InvClock_Id; i = helper::prev_index(i, R_Num))
				{
					r_queue.push_back(r_convex.convex_border_id[i]);
				}
				const uint16_t R = r_convex.convex_border_id[Tangent.R_InvClock_Id];
				r_queue.push_back(R);


				auto First = r_convex.convex_border_id.begin() + (Tangent.R_InvClock_Id + 1);
				auto Last = First + (r_convex.convex_border_id.size() - (Tangent.R_InvClock_Id + 1));
				r_convex.convex_border_id.erase(First, Last);
			}
			else
			{
				r_convex.right_is_line(Tangent.R_Clock_Id, Tangent.R_InvClock_Id, r_queue);
			}
			r_queue.shrink_to_fit();

			l_convex.r_border_id = l_convex.convex_border_id.size() + r_convex.r_border_id;
			vector_append(l_convex.convex_border_id, std::move(r_convex.convex_border_id));
			l_convex.is_line = false;
		}

	protected:
		struct convex_tangents
		{
			convex_tangents(const int32_t ClosetL, const int32_t ClosetR) //
				: R_Clock_Id(ClosetR)
				, R_InvClock_Id(ClosetR)
				, L_Clock_Id(ClosetL)
				, L_InvClock_Id(ClosetL)
			{}
			convex_tangents() //
				: R_Clock_Id(-1)
				, R_InvClock_Id(-1)
				, L_Clock_Id(-1)
				, L_InvClock_Id(-1)
			{}
			int32_t R_Clock_Id;
			int32_t R_InvClock_Id;
			int32_t L_Clock_Id;
			int32_t L_InvClock_Id;
		};

		template<typename T>
		static void vector_append(std::vector<T>& A, std::vector<T>&& B)
		{
			A.reserve(A.size() + B.size());
			A.insert(A.end(), std::make_move_iterator(B.begin()), std::make_move_iterator(B.end()));
			B.erase(B.begin(), B.end());
		}

		/* O(N^2)*/
		template<typename VectorGraph>
		static void find_closet_brute(
			const VectorGraph& graph,
			const convex_delaunay_2d<VectorType, IntType>& ConvexA,
			const convex_delaunay_2d<VectorType, IntType>& ConvexB,
			int32_t& ClosetA,
			int32_t& ClosetB)
		{

			double MinDistSquared = std::numeric_limits<double>::max();
			for (int32_t i = 0; i < ConvexA.convex_border_id.size(); ++i)
			{
				const VectorType& P = graph.Nodes[ConvexA.convex_border_id[i]];
				for (int32_t j = 0; j < ConvexB.convex_border_id.size(); j++)
				{
					const IntType NodeIdx = ConvexB.convex_border_id[j];
					const VectorType& Vec = P - graph.Nodes[NodeIdx];
					const double SizeSquared = (Vec[0] * Vec[0]) + (Vec[1] * Vec[1]); // 2d
					if (SizeSquared < MinDistSquared)
					{
						ClosetA = i;
						ClosetB = j;
						MinDistSquared = SizeSquared;
					}
				}
			}
		}

		static __forceinline bool is_C_Between_AB_cross0(const VectorType& A, const VectorType& B, const VectorType& C)
		{
			return C[0] >= std::min(A[0], B[0]) //
				&& C[0] <= std::max(A[0], B[0]) //
				&& C[1] >= std::min(A[1], B[1]) //
				&& C[1] <= std::max(A[1], B[1]);
		}

		template<typename VectorGraph>
		void find_not_clockwise_tangent(const VectorGraph& graph, const VectorType P, int32_t& BorderId) const
		{
			struct not_clockwise
			{
				__forceinline bool operator()(const VectorType& A, const VectorType& B, const VectorType& C) const
				{
					const double Cross = helper::cross2d<double, VectorType>(A, B, C);
					return !(Cross < 0.0) && (Cross > 0.0 || !is_C_Between_AB_cross0(A, C, B));
				};
			};

			int32_t PrevIdx = helper::prev_index(BorderId, convex_border_id.size());
			VectorType Vec0 = graph.Nodes[convex_border_id[BorderId]];
			VectorType VecPrev = graph.Nodes[convex_border_id[PrevIdx]];
			while (not_clockwise()(P, Vec0, VecPrev))
			{
				Vec0 = VecPrev;
				BorderId = PrevIdx;
				PrevIdx = helper::prev_index(PrevIdx, convex_border_id.size());
				VecPrev = graph.Nodes[convex_border_id[PrevIdx]];
			}
		}

		template<typename VectorGraph>
		void find_clockwise_tangent(const VectorGraph& VertPos, const VectorType VecIn, int32_t& BorderId) const
		{
			struct clockwise
			{
				__forceinline bool operator()(const VectorType& A, const VectorType& B, const VectorType& C) const
				{
					const double Cross = helper::cross2d<double, VectorType>(A, B, C);
					return !(Cross > 0.0) && (Cross < 0.0 || !is_C_Between_AB_cross0(A, C, B));
				};
			};

			int32_t Next = helper::next_index(BorderId, convex_border_id.size());
			VectorType Vec0 = VertPos.Nodes[convex_border_id[BorderId]];
			VectorType VecNext = VertPos.Nodes[convex_border_id[Next]];
			while (clockwise()(VecIn, Vec0, VecNext))
			{
				Vec0 = VecNext;
				BorderId = Next;
				Next = helper::next_index(Next, convex_border_id.size());
				VecNext = VertPos.Nodes[convex_border_id[Next]];
			}
		}

		/**
				 A_InvClock     B_Clock
					 |          |
				 _ _ o----------o_ _
				/     \        /     \
			   /       \      /       \
			  /         o    o         \
			  |         |    |    R    |
			  |    L    |    |         |
			  \         o    o         /
			   \       /      \       /
				\ _ _ o--------o _ _ /
					  |        |
				A_Clock       B_InvClock
		*/
		template<typename VectorGraph>
		static convex_tangents find_tangents_border_Id_to_convex(
			const VectorGraph& VertPos,
			const convex_delaunay_2d<VectorType, IntType>& ConvexA,
			const convex_delaunay_2d<VectorType, IntType>& ConvexB,
			int32_t ClosetA,
			int32_t ClosetB)
		{

			if (ClosetA == -1 || ClosetB == -1)
			{
				find_closet_brute(VertPos, ConvexA, ConvexB, ClosetA, ClosetB);
			}
			convex_tangents Tangent = convex_tangents(ClosetA, ClosetB);

			int32_t Old_A_InvClock;
			do
			{
				Old_A_InvClock = Tangent.L_InvClock_Id;
				VectorType V0 = VertPos.Nodes[ConvexA.convex_border_id[Tangent.L_InvClock_Id]];
				ConvexB.find_clockwise_tangent(VertPos, V0, Tangent.R_Clock_Id);

				V0 = VertPos.Nodes[ConvexB.convex_border_id[Tangent.R_Clock_Id]];
				ConvexA.find_not_clockwise_tangent(VertPos, V0, Tangent.L_InvClock_Id);
			} while (Old_A_InvClock != Tangent.L_InvClock_Id);

			int32_t Old_A_Clock;
			do
			{
				Old_A_Clock = Tangent.L_Clock_Id;
				VectorType V0 = VertPos.Nodes[ConvexA.convex_border_id[Tangent.L_Clock_Id]];
				ConvexB.find_not_clockwise_tangent(VertPos, V0, Tangent.R_InvClock_Id);

				V0 = VertPos.Nodes[ConvexB.convex_border_id[Tangent.R_InvClock_Id]];
				ConvexA.find_clockwise_tangent(VertPos, V0, Tangent.L_Clock_Id);
			} while (Old_A_Clock != Tangent.L_Clock_Id);

			return Tangent;
		}

		void left_is_line(const int32_t L_Clock_Id, const int32_t L_InvClock_Id, std::vector<IntType>& LQueue)
		{
			const int32_t L_Num = convex_border_id.size();
			if (L_Clock_Id == L_InvClock_Id)
			{
				LQueue.reserve(L_Num * 2 - 1);
				for (int32_t i = L_Num - 1; i > -L_Num; --i)
				{
					LQueue.push_back(*(convex_border_id.end() - 1 - std::abs(i)));
				}
				convex_border_id = {convex_border_id[L_Clock_Id]};
				l_border_id = 0;
			}
			else
			{
				if (L_Clock_Id != 0)
				{
					std::reverse(convex_border_id.begin(), convex_border_id.end());
					l_border_id = convex_border_id.size() - 1;
				}
				else
				{
					l_border_id = 0;
				}

				for (int32_t i = L_Num - 1; i >= 0; --i)
				{
					LQueue.push_back(convex_border_id[i]);
				}
			}
		}
		void right_is_line(const int32_t R_Clock_Id, const int32_t R_InvClock_Id, std::vector<IntType>& RQueue)
		{
			const int32_t R_Num = convex_border_id.size();
			if (R_Clock_Id == R_InvClock_Id)
			{
				RQueue.reserve(R_Num * 2 - 1);
				for (int32_t i = R_Num - 1; i > -R_Num; --i)
				{
					RQueue.push_back(convex_border_id[std::abs(i)]);
				}
				convex_border_id = {convex_border_id[R_Clock_Id]};
				r_border_id = 0;
			}
			else
			{
				if (R_Clock_Id != 0)
				{
					std::reverse(convex_border_id.begin(), convex_border_id.end());
					r_border_id = 0;
				}
				else
				{
					r_border_id = convex_border_id.size() - 1;
				}

				for (int32_t i = 0; i < R_Num; ++i)
				{
					RQueue.push_back(convex_border_id[i]);
				}
			}
		}
	};


	template<typename VectorType, typename IntType = uint16_t, typename LinkType = std::vector<IntType>>
	class vector_graph_triangulator
	{
	public:
		using GraphType = vector_graph<VectorType, IntType, LinkType>;

		GraphType& Graph;

		vector_graph_triangulator(GraphType& InGraph) : Graph(InGraph)
		{
			//triangulate_2d();
		}


		__forceinline void triangulate_2d() { triangulate2d_impl(); }

	private:
		bool sort_adjacent(const IntType InNodeID, const bool bSortNoBorder = false)
		{
			const VectorType& V = Graph.Nodes[InNodeID];
			LinkType& CurrentLinks = Graph.NodeLinks[InNodeID];

			if (CurrentLinks.size() < 2)
			{
				return true;
			}
			if (CurrentLinks.size() == 2) //todo
			{
				if (!(helper::cross2d<double, VectorType>(V, Graph.Nodes[CurrentLinks[0]], Graph.Nodes[CurrentLinks[1]]) >= 0.0))
				{
					std::swap(CurrentLinks[0], CurrentLinks[1]);
				}
				return true;
			}

			std::vector<std::pair<IntType, float>> Angles;
			Angles.reserve(CurrentLinks.size());
			for (const auto& It : CurrentLinks)
			{
				const VectorType V1 = Graph.Nodes[It] - V;
				float A = std::atan2(V1[0], V1[1]);
				while (A > M_PI)
				{
					A -= (M_PI * 2.0f);
				}
				while (A < -M_PI)
				{
					A += (M_PI * 2.0f);
				}
				Angles.push_back(std::pair<IntType, float>(It, A));
			}
			std::sort(
				Angles.begin(),
				Angles.end(),
				[](const std::pair<IntType, float>& A, const std::pair<IntType, float>& B) //
				{																		   //
					return A.second > B.second;
				});


			float Max = std::numeric_limits<float>::min();
			int32_t MaxId = -1;
			int32_t CloseId = 0;
			for (int32_t i = 0; i < Angles.size(); ++i)
			{
				CurrentLinks[i] = Angles[i].first;
				const int32_t Next = helper::next_index(i, Angles.size());

				const float f = Next != 0 //
					? Angles[i].second - Angles[Next].second
					: Angles.back().second - Angles[Next].second + (2.0 * M_PI);

				if (f > Max)
				{
					Max = f;
					MaxId = Next;
				}
				if (Angles[i].first > InNodeID && Angles[i].first < Angles[CloseId].first)
				{
					CloseId = i;
				}
			}

			if (Max >= M_PI)
			{
				helper::rotate_to_index(CurrentLinks, MaxId);
				return true;
			}
			if (bSortNoBorder)
			{
				helper::rotate_to_index(CurrentLinks, CloseId);
			}
			return false;
		}

		std::vector<IntType> sorted_2d_id() const
		{
			std::vector<IntType> Out;

			const int32_t Num = Graph.Nodes.size();
			Out.reserve(Num);
			for (int32_t i = 0; i < Num; ++i)
			{
				Out.push_back(static_cast<IntType>(i));
			}

			std::sort(
				Out.begin(),
				Out.end(),
				[&](const IntType A, const IntType B)
				{
					const VectorType& V1 = Graph.Nodes[A];
					const VectorType& V2 = Graph.Nodes[B];
					return V1[0] < V2[0] || (V1[0] == V2[0] && V1[1] < V2[1]);
				});

			return std::move(Out);
		}

		convex_delaunay_2d<VectorType, IntType> triangle_small(const helper::range_view<IntType> Px)
		{
			Graph.add_edge(Px[0], Px[1]);

			convex_delaunay_2d<VectorType, IntType> Out;
			Out.l_border_id = 0;
			Out.r_border_id = 1;
			Out.convex_border_id.reserve(Px.size());
			Out.convex_border_id.push_back(Px[0]);
			Out.convex_border_id.push_back(Px[1]);

			if (Px.size() == 3)
			{
				Out.convex_border_id.push_back(Px[2]);
				Graph.add_edge(Px[1], Px[2]);

				const double Cross = helper::cross2d<double, VectorType>(Graph.Nodes[Px[0]], Graph.Nodes[Px[1]], Graph.Nodes[Px[2]]);
				if (Cross == 0.0)
				{
					Out.is_line = true;
				}
				else
				{
					Graph.add_edge(Px[0], Px[2]);
				}

				if (Cross >= 0.0)
				{
					Out.r_border_id = 2;
				}
				else
				{
					std::swap(Out.convex_border_id[1], Out.convex_border_id[2]);
				}
			}
			else if (Px.size() == 2)
			{
				Out.is_line = true;
			}
			return Out;
		}


		/** Delaunay triangulation */
		void delaunay_corridor_2d(std::vector<IntType>&& LQueue, std::vector<IntType>&& RQueue, bool bConst_L = false, bool bConst_R = false)
		{
			/**
				  LArr[0]----->RArr[0]
						\      /
					LArr[i]   RArr[j]
						 |    |
						 |    |
					  LNext  RNext
						/      \
			 LArr.Last()<-------RArr.Last()
			*/

			IntType L = LQueue.back();
			LQueue.pop_back();
			IntType R = RQueue.back();
			RQueue.pop_back();

			Graph.add_edge(L, R);

			VectorType LvNext, RvNext;

			using FCircle_Rsq = helper::circle2d_rsq<double, VectorType>;
			FCircle_Rsq LCircle, RCircle;

			while (LQueue.size() > 0 || RQueue.size() > 0)
			{
				IntType L_Next = LQueue.size() > 0 ? LQueue.back() : std::numeric_limits<IntType>::max();
				IntType R_Next = RQueue.size() > 0 ? RQueue.back() : std::numeric_limits<IntType>::max();
				const VectorType& Lv = Graph.Nodes[L];
				const VectorType& Rv = Graph.Nodes[R];
				bool bLeft = false;
				bool bRight = false;

				if (L_Next != std::numeric_limits<IntType>::max())
				{
					LvNext = Graph.Nodes[L_Next];
					bLeft = helper::cross2d<double, VectorType>(Rv, Lv, LvNext) > 0.0;
					if (bLeft)
					{
						LCircle = FCircle_Rsq(Lv, Rv, LvNext);
						if (!bConst_L) //todo delaunay_corridor_2d bConst_L
						{
							IntType OppID;
							VectorType OppV;
							do
							{
								OppID = std::numeric_limits<IntType>::max();
								const LinkType& NextLinks = Graph.NodeLinks[L_Next];
								const LinkType& LLin = Graph.NodeLinks[L];
								for (const IntType Link : LLin)
								{
									if (Link != R && GraphType::contains_link(NextLinks, Link))
									{
										const VectorType& LinkV = Graph.Nodes[Link];
										if (helper::cross2d<double, VectorType>(LinkV, LvNext, Lv) > 0.0)
										{
											if (OppID == std::numeric_limits<IntType>::max()				   //check is inside triangle
												|| (helper::cross2d<double, VectorType>(LinkV, Lv, OppV) > 0.0 //
													&& helper::cross2d<double, VectorType>(LinkV, OppV, LvNext) > 0.0))
											{
												OppV = LinkV;
												OppID = Link;
											}
										}
									}
								}
								if (OppID != std::numeric_limits<IntType>::max() //
									&& LCircle.IsInCircle(OppV)					 //
									&& helper::cross2d<double, VectorType>(Rv, Lv, OppV) > 0.0)
								{
									Graph.remove_edge(L, L_Next);
									LQueue.push_back(OppID);
									L_Next = OppID;
									LvNext = OppV;
									LCircle = FCircle_Rsq(Lv, Rv, LvNext);
								}
								else
								{
									OppID = std::numeric_limits<IntType>::max();
								}
							} while (OppID != std::numeric_limits<IntType>::max());
						}
					}
				}
				if (R_Next != std::numeric_limits<IntType>::max())
				{
					RvNext = Graph.Nodes[R_Next];
					bRight = helper::cross2d<double, VectorType>(Lv, Rv, RvNext) < 0.0;
					if (bRight)
					{
						RCircle = FCircle_Rsq(Lv, Rv, RvNext);
						if (!bConst_R) //todo delaunay_corridor_2d bConst_R
						{
							IntType OppID;
							VectorType OppV;
							do
							{
								OppID = std::numeric_limits<IntType>::max();
								const LinkType& NextLinks = Graph.NodeLinks[R_Next];
								const LinkType& RLin = Graph.NodeLinks[R];
								for (const IntType Link : RLin)
								{
									if (Link != L && GraphType::contains_link(NextLinks, Link))
									{
										const VectorType& LinkV = Graph.Nodes[Link];
										if (helper::cross2d<double, VectorType>(LinkV, Rv, RvNext) > 0.0)
										{
											if (OppID == std::numeric_limits<IntType>::max()					   //check is inside triangle
												|| (helper::cross2d<double, VectorType>(LinkV, RvNext, OppV) > 0.0 //
													&& helper::cross2d<double, VectorType>(LinkV, OppV, Rv) > 0.0))
											{
												OppV = LinkV;
												OppID = Link;
											}
										}
									}
								}
								if (OppID != std::numeric_limits<IntType>::max() //
									&& RCircle.IsInCircle(OppV)					 //
									&& helper::cross2d<double, VectorType>(Lv, Rv, OppV) < 0.0)
								{
									Graph.remove_edge(R, R_Next);
									RQueue.push_back(OppID);
									R_Next = OppID;
									RvNext = OppV;
									RCircle = FCircle_Rsq(Lv, Rv, RvNext);
								}
								else
								{
									OppID = std::numeric_limits<IntType>::max();
								}
							} while (OppID != std::numeric_limits<IntType>::max());
						}
					}
				}

				//if (!bLeft && !bRight) return;

				if (bLeft && bRight)
				{
					if (LCircle.IsInCircle(RvNext))
					{
						bLeft = false;
					}
					else if (RCircle.IsInCircle(LvNext))
					{
						bRight = false;
					}
					else
					{
						LCircle.RadiusSquared < RCircle.RadiusSquared ? (bRight = false) : (bLeft = false);
					}
				}

				if (bLeft)
				{
					Graph.add_edge(L_Next, R);
					LQueue.pop_back();
					L = L_Next;
				}
				else
				{
					Graph.add_edge(R_Next, L);
					RQueue.pop_back();
					R = R_Next;
				}
			}
		}

		__forceinline void union_convex(convex_delaunay_2d<VectorType, IntType>& LConvex, convex_delaunay_2d<VectorType, IntType>&& RConvex)
		{
			std::vector<IntType> LQueue, RQueue;
			LConvex.convex_union(std::move(RConvex), Graph, LQueue, RQueue);
			delaunay_corridor_2d(std::move(LQueue), std::move(RQueue));
		}

		std::vector<convex_delaunay_2d<VectorType, IntType>> triangulator_2d_small()
		{
			std::vector<convex_delaunay_2d<VectorType, IntType>> Out;

			std::vector<IntType> SortId = sorted_2d_id();
			helper::range_view<IntType> Px(&SortId[0], SortId.size());

			const int32_t Last = Px.size() % 3;
			const int32_t Num = Px.size() / 3 + (Last > 0 ? 1 : 0);
			Out.resize(Num, convex_delaunay_2d<VectorType, IntType>());

			const int32_t LoopNum = Px.size() / 3 - (Last == 1 ? 1 : 0);

			concurrency::parallel_for(
				0,
				LoopNum,
				1,
				[&](const int32_t Idx) //
				{					   //
					Out[Idx] = triangle_small(Px.slice(Idx * 3, 3));
				});

			if (Last > 0)
			{
				if (Last == 1)
				{
					Out[LoopNum] = triangle_small(Px.slice((Px.size() - 4), 2));
					Out[LoopNum + 1] = triangle_small(Px.slice((Px.size() - 2), 2));
				}
				else
				{
					Out[LoopNum] = triangle_small(Px.slice((Px.size() - Last), Last));
				}
			}
			return std::move(Out);
		}

		void triangulate2d_impl()
		{
			std::vector<convex_delaunay_2d<VectorType, IntType>> ConvexArray = triangulator_2d_small();
			int32_t Num = ConvexArray.size();

			while (Num > 1)
			{
				const int32_t OneLast = Num % 2;
				Num = Num / 2;

				concurrency::parallel_for(
					0,
					Num,
					1,
					[&](int32_t Idx)
					{
						Idx *= 2;
						union_convex(ConvexArray[Idx], std::move(ConvexArray[Idx + 1]));
					});

				Num += OneLast;
				for (int32_t i = 1; i < Num; ++i)
				{
					std::swap(ConvexArray[i], ConvexArray[i * 2]);
				}
				//ConvexArray.resize(Num); // free memory
			}
			ConvexArray.erase(ConvexArray.begin(), ConvexArray.end());

			// for voronoi graph
			concurrency::parallel_for(
			0,
			this->Nodes.size(),
			1,
			[&](const int32_t Idx) //
			{					   //
				this->sort_adjacent(static_cast<IntType>(Idx), false);
			});

		}
	};


	template< //
		typename VectorType,
		typename IntType = uint16_t,
		typename LinkType = std::vector<IntType>>
	static vector_graph<VectorType, IntType, LinkType> triangulate(VectorType* VecBegin, int32_t Num)
	{
		using GraphType = vector_graph<VectorType, IntType, LinkType>;
		GraphType Graph(VecBegin, Num);
		auto Triangulator = vector_graph_triangulator<VectorType, IntType, LinkType>(Graph);
		Triangulator.triangulate_2d();
		return std::move(Graph);
	}

	template< //
		typename ArrType,
		typename VectorType,
		typename IntType = uint16_t,
		typename LinkType = std::vector<IntType>>
	static vector_graph<VectorType, IntType, LinkType> triangulate(ArrType& In)
	{
		auto IterBegin = In.begin();
		int32_t Size = std::distance(IterBegin, In.end());
		return triangulate<VectorType, IntType, LinkType>(IterBegin, Size);
	}

	template< //
		typename VectorType,
		typename IntType = uint16_t,
		typename LinkType = std::vector<IntType>>
	static vector_graph<VectorType, IntType, LinkType> triangulate(std::vector<VectorType>& In)
	{
		return triangulate<std::vector<VectorType>, VectorType, IntType, LinkType>(In);
	}


	template< //
		typename VectorType,
		typename IntType = uint16_t,
		typename LinkType = std::vector<IntType>>
	static void triangulate_graph(vector_graph<VectorType, IntType, LinkType>& Graph)
	{
		vector_graph_triangulator<VectorType, IntType, LinkType>(Graph).triangulate_2d();
	}

} // namespace Triangulate