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

#ifndef _4ti2_qsolve__PlusTree_
#define _4ti2_qsolve__PlusTree_

#include <vector>

namespace _4ti2_
{

template <class IndexSet>
class PlusTree
{
public:
    PlusTree();
    ~PlusTree();

    void insert(const std::vector<IndexSet>& supports);
    bool dominated(const IndexSet& b, int index1, int index2) const;
    void dump() const;
    void find_diff(std::vector<int>& indices, const IndexSet& supp, int diff) const;
    void insert(const IndexSet& support, int index);

    void clear();

private:
    class TreeNode
    {
    public:
        TreeNode() {}
        virtual ~TreeNode() {}
        virtual TreeNode* insert(const IndexSet& support, int start, int index) = 0;
        virtual bool dominated(const IndexSet& b, int index1, int index2) const = 0;
        virtual void find_diff(std::vector<int>& indices, const IndexSet& supp, 
                                int diff, int orig_diff) const = 0;
        virtual void dump(Index level, Index index) const = 0;
    };

    class TreeLeaf : public TreeNode {
    public:
        TreeLeaf(Index _index, const IndexSet& _is) : i(_index), is(_is) {}
        ~TreeLeaf() {}
        TreeNode* insert(const IndexSet& support, int start, int index);
        bool dominated(const IndexSet& b, int index1, int index2) const {
            if (is.set_subset(b) && i != index1 && i != index2) { return true; }
            return false;
        }
        void find_diff(std::vector<int>& indices, const IndexSet& supp, int diff, int orig_diff) const {
            IndexSet tmp(is);
            tmp.set_intersection(supp);
            if (tmp.singleton()) { indices.push_back(i); }
        }
        void dump(Index level, Index index) const {
            std::string indent(level*INDENT, ' ');
            *out << indent << "B[" << i << "]: ";
            for (Index i = index; i < is.get_size(); ++i) {
                if (is[i]) { *out << i << " "; }
            }
            *out << "\n";
        }
    private:
        Index i;
        IndexSet is;
    };

    class TreeBranch : public TreeNode {
    public:
        TreeBranch();
        ~TreeBranch();
        TreeNode* insert(const IndexSet& support, int start, int index);
        bool dominated(const IndexSet& b, int index1, int index2) const;
        void find_diff(std::vector<int>& indices, const IndexSet& supp, int diff, int orig_diff) const;
        void dump(Index level, Index index) const {
            std::string indent(level*INDENT, ' ');
            *out << indent << "Limb:\n";
            for (unsigned i = 0; i < nodes.size(); ++i) {
                *out << indent << "-" << nodes[i].first << "\n";
                nodes[i].second->dump(level+1, nodes[i].first+1);
            }
        }
        void add(Index index, TreeNode* node) { nodes.push_back(std::pair<int,TreeNode*>(index,node)); }
        void clear();
    private:
        std::vector<std::pair<int,TreeNode*> > nodes;
    };

    static const int INDENT;

    TreeBranch root;
};

template <class IndexSet>
PlusTree<IndexSet>::PlusTree()
{
}

template <class IndexSet>
PlusTree<IndexSet>::~PlusTree()
{
}

template <class IndexSet>
inline
void
PlusTree<IndexSet>::insert(const std::vector<IndexSet>& supports)
{
    for (unsigned int i = 0; i < supports.size(); ++i) {
        insert(supports[i], i);
    }
}

template <class IndexSet>
inline
void
PlusTree<IndexSet>::insert(const IndexSet& support, int index)
{
    root.insert(support, 0, index);
}

template <class IndexSet>
inline
bool
PlusTree<IndexSet>::dominated(const IndexSet& support, int index1, int index2) const
{
    return root.dominated(support, index1, index2);
}

template <class IndexSet>
inline
void
PlusTree<IndexSet>::dump() const
{
    *out << "Plus Tree Dump:\n";
    root.dump(0,0);
}

template <class IndexSet>
inline
void
PlusTree<IndexSet>::clear()
{
    root.clear();
}

template <class IndexSet>
inline
void
PlusTree<IndexSet>::find_diff(std::vector<int>& indices, const IndexSet& supp, int diff) const
{
    root.find_diff(indices, supp, diff, diff);
}

} // namespace _4ti2_

// Include template defitinions.
#include "qsolve/PlusTree.hpp"

#endif
