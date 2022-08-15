// Copyright 2022 Alexandr Marchenko. All Rights Reserved.

#pragma once

#include <vector>
#include <algorithm>

#define _USE_MATH_DEFINES
#include <cmath>


namespace helper
{
	template<typename Real = double, typename VectorType>
	static __forceinline Real cross2d(const VectorType& A, const VectorType& B, const VectorType& C)
	{
		return																											//
			(static_cast<Real>(B[0]) - static_cast<Real>(A[0])) * (static_cast<Real>(C[1]) - static_cast<Real>(A[1])) - //
			(static_cast<Real>(B[1]) - static_cast<Real>(A[1])) * (static_cast<Real>(C[0]) - static_cast<Real>(A[0]));
	}

	template<typename ElementType>
	struct range_view
	{
		ElementType* m_data_ptr;
		int32_t m_num;
		range_view<ElementType>() : m_data_ptr(nullptr), m_num(0) {}
		range_view<ElementType>(ElementType* Ptr, int32_t InNum) : m_data_ptr(Ptr), m_num(InNum) {}

		range_view<ElementType> slice(const int32_t Index, const int32_t InNum) const { return range_view<ElementType>(m_data_ptr + Index, InNum); }
		int32_t size() const { return m_num; }
		ElementType& operator[](const uint32_t Index) const { return m_data_ptr[Index]; }

		ElementType* begin() const { return m_data_ptr; }
		ElementType* end() const { return m_data_ptr + m_num; }
	};

} // namespace helper


/** Oriented Vector Graph*/
template< //
	typename VectorType,
	typename IntType = uint16_t,
	typename LinkType = std::vector<IntType>>
class vector_graph
{
public:
	using TriangleType = std::tuple<IntType, IntType, IntType>;
	using EdgeType = std::tuple<IntType, IntType>;
	using SelfType = vector_graph<VectorType, IntType, LinkType>;

	/*****************************/

	/** range Nodes */
	helper::range_view<VectorType> Nodes;
	/** Nodes Links */
	std::vector<LinkType> NodeLinks;

	/*****************************/

	vector_graph<VectorType, IntType, LinkType>() {}
	vector_graph<VectorType, IntType, LinkType>(VectorType* VecBegin, int32_t Num) : Nodes(VecBegin, Num)
	{
		if (Nodes.size() > 0)
		{
			this->NodeLinks.resize(this->Nodes.size(), LinkType());
		}
	}
	vector_graph<VectorType, IntType, LinkType>(std::vector<VectorType>& In) : Nodes(&In[0], In.size())
	{
		if (Nodes.size() > 0)
		{
			this->NodeLinks.resize(this->Nodes.size(), LinkType());
		}
	}
	~vector_graph<VectorType, IntType, LinkType>() = default;

	/*****************************/

	//todo add_unique check
	void add_edge(const IntType From, const IntType To, const bool bTwoSided = true)
	{
		NodeLinks[From].push_back(To);
		if (bTwoSided)
		{
			NodeLinks[To].push_back(From);
		}
	}
	void remove_edge(const IntType U, const IntType V)
	{
		auto& ConU = NodeLinks[U];
		auto It1 = std::find(ConU.begin(), ConU.end(), V);
		if (It1 != ConU.end())
		{
			ConU.erase(It1);
		}

		auto& ConV = NodeLinks[V];
		auto It2 = std::find(ConV.begin(), ConV.end(), U);
		if (It1 != ConV.end())
		{
			ConV.erase(It2);
		}
	}

	/*****************************/


	static __forceinline bool contains_link(const LinkType& InLinks, const IntType L) { return InLinks.end() != std::find(InLinks.begin(), InLinks.end(), L); }
	bool contains_link(const IntType ContainerIdx, const IntType L) const { return SelfType::contains_link(NodeLinks[ContainerIdx], L); }

	bool is_connected(const IntType A, const IntType B, const bool bTwoSided = false) const
	{
		if (!bTwoSided)
		{
			return (NodeLinks[A].size() < NodeLinks[B].size()) //
				? contains_link(A, B)						   //
				: contains_link(B, A);						   //
		}
		return contains_link(A, B) || contains_link(B, A);
	}
	bool is_triangle(const IntType A, const IntType B, const IntType C, const bool bTwoSided = false) const
	{
		return is_connected(A, B, bTwoSided) && is_connected(A, C, bTwoSided) && is_connected(B, C, bTwoSided);
	}


	template<typename LambdaType>
	void for_each_edge(LambdaType Lambda)
	{
		for (int32_t i = 0; i < Nodes.size(); ++i)
		{
			VectorType& V1 = Nodes[i];
			for (auto j : NodeLinks[i])
			{
				VectorType& V2 = Nodes[j];
				if (V1[0] < V2[0] || (V1[0] == V2[0] && V1[1] < V2[1]))
				{
					Lambda(V1, V2);
				}
			}
		}
	};
	template<typename LambdaType>
	void for_each_edge(LambdaType Lambda) const
	{
		for (int32_t i = 0; i < Nodes.size(); ++i)
		{
			const VectorType& V1 = Nodes[i];
			for (auto j : NodeLinks[i])
			{
				const VectorType& V2 = Nodes[j];
				if (V1[0] < V2[0] || (V1[0] == V2[0] && V1[1] < V2[1]))
				{
					Lambda(V1, V2);
				}
			}
		}
	};

	// todo: works only if have sorted links
	bool IsBorder(IntType InNodeID) const
	{
		const LinkType& CurrentLinks = NodeLinks[InNodeID];

		if (CurrentLinks.size() == 2) return true;

		const double Cross = helper::cross2d<double, VectorType>(Nodes[InNodeID], Nodes[CurrentLinks.back()], Nodes[CurrentLinks.front()]);
		return Cross <= 0.0;
	}

	std::vector<TriangleType> get_connected_triangles(IntType P) const //todo get_connected_triangles
	{
		std::vector<TriangleType> Out;
		const auto& Links = NodeLinks[P];
		const int32_t Size = Links.size();
		Out.reserve(Size);
		if (Size > 2)
		{
			for (int32_t i = 0; i < Size; ++i)
			{
				for (int32_t j = 1; j < Size; ++j)
				{
					if (is_connected(Links[i], Links[j]))
					{
						Out.push_back(std::make_tuple(P, Links[i], Links[j]));
					}
				}
			}
		}
		return Out;
	}
	std::vector<EdgeType> get_connected_edges(IntType P) const //todo get_connected_edges
	{
		std::vector<EdgeType> Out;
		const auto& Links = NodeLinks[P];
		const int32_t Size = Links.size();
		Out.reserve(Size);
		if (Size > 2)
		{
			for (int32_t i = 0; i < Size; i++)
			{
				for (int32_t j = 1; j < Size; j++)
				{
					if (is_connected(Links[i], Links[j]))
					{
						Out.push_back(std::make_tuple(Links[i], Links[j]));
					}
				}
			}
		}
		return Out;
	}
};
