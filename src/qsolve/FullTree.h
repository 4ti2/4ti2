/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Peter Malkin.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 
*/

#ifndef _4ti2_qsolve__FullTree_
#define _4ti2_qsolve__FullTree_

#include <vector>

namespace _4ti2_
{

template <class IndexSet>
class FullTree
{
public:
    FullTree();
    ~FullTree();

    void insert(const IndexSet& support, Index index);
    void insert(const std::vector<IndexSet>& supports);
    void insert(const std::vector<IndexSet>& supports, Index start, Index end);
    bool dominated(const IndexSet& b, Index index1, Index index2) const;
    void dump() const;
    void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const;
    void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const;
    void print_statistics();

    void clear();

private:
    static long long int num_nodes;
    static long long int num_limbs;
    static long long int num_buds;
    static long long int num_buckets;
    static const int INDENT;

    class TreeNode
    {
    public:
        TreeNode() {}
        virtual ~TreeNode() {}
        virtual TreeNode* insert(const IndexSet& supp, const IndexSet& parent, Index index) = 0;
        virtual bool dominated(const IndexSet& b, Index index1, Index index2) const = 0;
        virtual void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const = 0;
        virtual void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const = 0;
        virtual void dump(Index level) const = 0;

        virtual void compute_statistics() = 0;
    protected:
    };

    class TreeBud: public TreeNode {
    public:
        TreeBud(Index _index) : i(_index) {}
        virtual ~TreeBud() {}
        virtual TreeNode* insert(const IndexSet& supp, const IndexSet& parent, Index index);
        virtual bool dominated(const IndexSet& b, Index index1, Index index2) const;
        virtual void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const;
        virtual void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const;
        virtual void dump(Index level) const;
        virtual void compute_statistics() { ++num_nodes; ++num_buds; }
    private:
        Index i;
    };

    class TreeLimb : public TreeNode {
    public:
        TreeLimb();
        virtual ~TreeLimb();
        virtual TreeNode* insert(const IndexSet& supp, const IndexSet& parent, Index index);
        virtual void insert(const IndexSet& supp, Index index, TreeNode* node);
        virtual bool dominated(const IndexSet& b, Index index1, Index index2) const;
        virtual void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const;
        virtual void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const;
        virtual void clear();

        virtual void dump(Index level) const;
        virtual void compute_statistics();
    private:
        std::vector<IndexSet> supps;
        std::vector<TreeNode*> nodes;
        std::vector<Index> indices;
    };

    class TreeBucket: public TreeNode {
    public:
        TreeBucket();
        ~TreeBucket();
        virtual TreeNode* insert(const IndexSet& supp, const IndexSet& parent, Index index);
        virtual bool dominated(const IndexSet& b, Index index1, Index index2) const;
        virtual void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const;
        virtual void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const;
        virtual void dump(int level) const;
        virtual void compute_statistics() { ++num_nodes; ++num_buckets; }
    private:
        std::vector<IndexSet> supps;
        std::vector<Index> indices;
        static const Size MAX_BUCKET_SIZE;
    };

    TreeLimb root;
};

} // namespace _4ti2_

// Include template defitinions.
#include "qsolve/FullTree.hpp"

#endif
