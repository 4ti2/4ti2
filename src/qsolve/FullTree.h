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

    void insert(const std::vector<IndexSet>& supports);
    void insert(const std::vector<IndexSet>& supports, Index start, Index end);
    bool dominated(const IndexSet& b, Index index1, Index index2) const;
    void dump() const;
    void find_diff(std::vector<Index>& inds, const IndexSet& supp, Size diff) const;
    void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const;
    void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const;
    void find(std::vector<Index>& inds, std::vector<Index>& inds0, const IndexSet& zeros, const IndexSet& supp, Size count) const;
    void insert(const IndexSet& support, Index index);
    void print_statistics();

    void clear();

private:
    class TreeNode
    {
    public:
        TreeNode() {}
        virtual ~TreeNode() {}
        virtual TreeNode* insert(const IndexSet& support, Index start, Index index) = 0;
        virtual bool dominated(const IndexSet& b, Index index1, Index index2) const = 0;
        virtual void find_diff(std::vector<Index>& inds, const IndexSet& supp, Size diff, Size orig_diff) const = 0;
        virtual void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const = 0;
        virtual void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const = 0;
        virtual void find(std::vector<Index>& inds, std::vector<Index>& inds0, const IndexSet& zeros, const IndexSet& supp, Size count) const = 0;
        virtual void dump(Index level, Index index) const = 0;
        //virtual TreeNode* finalise() = 0;

        virtual void compute_statistics() = 0;
        static void print_statistics();
        static void reset_statistics();
    protected:
        static long long int num_nodes;
        static long long int num_branches;
        static long long int num_limbs;
        static long long int num_arms;
        static long long int num_leaves;
        static long long int num_buds;
        static long long int num_singletons;
    };

#if 0
    class TreeLeaf : public TreeNode {
    public:
        TreeLeaf(Index _index, const IndexSet& _is) : i(_index), is(_is) {}
        virtual ~TreeLeaf() {}
        virtual TreeNode* insert(const IndexSet& support, Index start, Index index);
        virtual bool dominated(const IndexSet& b, Index index1, Index index2) const;
        virtual void find_diff(std::vector<Index>& inds, const IndexSet& supp, Index diff, Index orig_diff) const;
        virtual void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const;
        virtual void dump(Index level, Index index) const;
        virtual void compute_statistics() { ++TreeNode::num_nodes; ++TreeNode::num_leaves; }
    private: 
        Index i;
        IndexSet is;
    };
#endif

    class TreeBud: public TreeNode {
    public:
        TreeBud(Index _index, const IndexSet& _is) : i(_index), is(_is) {}
        virtual ~TreeBud() {}
        virtual TreeNode* insert(const IndexSet& support, Index start, Index index);
        virtual bool dominated(const IndexSet& b, Index index1, Index index2) const;
        virtual void find_diff(std::vector<Index>& inds, const IndexSet& supp, Index diff, Index orig_diff) const;
        virtual void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const;
        virtual void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const;
        virtual void find(std::vector<Index>& inds, std::vector<Index>& inds0, const IndexSet& zeros, const IndexSet& supp, Size count) const;
        virtual void dump(Index level, Index index) const;
        virtual void compute_statistics() { ++TreeNode::num_nodes; ++TreeNode::num_buds; }
    private:
        Index i;
        IndexSet is;
    };

#if 0
    class TreeBranch : public TreeNode {
    public:
        TreeBranch();
        TreeBranch(Index next, TreeNode* node);
        virtual ~TreeBranch();
        virtual TreeNode* insert(const IndexSet& support, Index start, Index index);
        virtual bool dominated(const IndexSet& b, Index index1, Index index2) const;
        virtual void find_diff(std::vector<Index>& inds, const IndexSet& supp, Index diff, Index orig_diff) const;
        virtual void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const;
        virtual void clear();

        virtual void dump(Index level, Index index) const;
        virtual void compute_statistics();
    private:
        std::vector<std::pair<Index,TreeNode*> > nodes;
    };
#endif

    class TreeLimb : public TreeNode {
    public:
        TreeLimb(Index start, const IndexSet& s1, TreeNode* n1, const IndexSet& s2, TreeNode* n2);
        virtual ~TreeLimb();
        virtual TreeNode* insert(const IndexSet& support, Index start, Index index);
        virtual bool dominated(const IndexSet& b, Index index1, Index index2) const;
        virtual void find_diff(std::vector<Index>& inds, const IndexSet& supp, Index diff, Index orig_diff) const;
        virtual void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const;
        virtual void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const;
        virtual void find(std::vector<Index>& inds, std::vector<Index>& inds0, const IndexSet& zeros, const IndexSet& supp, Size count) const;
        virtual void clear();

        virtual void dump(Index level, Index index) const;
        virtual void compute_statistics();
    private:
        std::vector<Index> indices;
        std::vector<IndexSet> supps;
        std::vector<TreeNode*> nodes;
        IndexSet supp_intersection;
    };

    class TreeArm : public TreeNode {
    public:
        TreeArm();
        TreeArm(const IndexSet& support, Index next, TreeNode* node);
        TreeArm(Index start, const IndexSet& s1, TreeNode* n1, const IndexSet& s2, TreeNode* n2);
        virtual ~TreeArm();
        virtual TreeNode* insert(const IndexSet& support, Index start, Index index);
        virtual bool dominated(const IndexSet& b, Index index1, Index index2) const;
        virtual void find_diff(std::vector<Index>& inds, const IndexSet& supp, Index diff, Index orig_diff) const;
        virtual void find_singleton_diff(std::vector<Index>& inds, const IndexSet& supp) const;
        virtual void find(std::vector<Index>& inds, const IndexSet& zeros, const IndexSet& supp, Size count) const;
        virtual void find(std::vector<Index>& inds, std::vector<Index>& inds0, const IndexSet& zeros, const IndexSet& supp, Size count) const;
        virtual void clear();

        virtual void dump(Index level, Index index) const;
        virtual void compute_statistics();
    private:
        std::vector<Index> indices;
        std::vector<IndexSet> supps;
        std::vector<TreeNode*> nodes;
    };

    static const int INDENT;

    //TreeBranch root;
    //TreeLimb root;
    TreeArm root;
};

} // namespace _4ti2_

// Include template defitinions.
#include "qsolve/FullTree.hpp"

#endif
