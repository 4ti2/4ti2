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

#ifndef _4ti2_groebner__OnesTree_
#define _4ti2_groebner__OnesTree_

#include <vector>

namespace _4ti2_
{

template <class IndexSet>
class OnesTreeNode
{
public:
    OnesTreeNode();
    virtual ~OnesTreeNode();

    virtual void insert(const IndexSet& support, int start, int remaining, int index) = 0;
    virtual void dump(int level) = 0;
    virtual bool dominated(const IndexSet& b, int index1, int index2) = 0;
    virtual int index_dominated(const IndexSet& b, int index1, int index2) = 0;
    virtual int find(const IndexSet& b, int start, int remaining) = 0;
    virtual void find_diff(std::vector<int>& indices, const IndexSet& supp, int diff) = 0;
    virtual void find_diff(std::vector<int>& indices, const IndexSet& supp1, int diff1, const IndexSet& supp2, int diff2) = 0;

protected:
    static const int INDENT;
};

template <class IndexSet>
class OnesTreeBranch : public OnesTreeNode<IndexSet>
{
public:
    OnesTreeBranch();
    virtual ~OnesTreeBranch();

    virtual void insert(const IndexSet& support, int start, int remaining, int index);
    virtual void dump(int level);
    virtual bool dominated(const IndexSet& b, int index1, int index2);
    virtual int index_dominated(const IndexSet& b, int index1, int index2);
    virtual int find(const IndexSet& b, int start, int remaining);
    virtual void find_diff(std::vector<int>& indices, const IndexSet& supp, int diff);
    virtual void find_diff(std::vector<int>& indices, const IndexSet& supp1, int diff1, const IndexSet& supp2, int diff2);

protected:
    std::vector<std::pair<int,OnesTreeNode<IndexSet>*> >nodes;
};

template <class IndexSet>
class OnesTreeLeaf : public OnesTreeNode<IndexSet>
{
public:
    OnesTreeLeaf();
    virtual ~OnesTreeLeaf();
 
    virtual void insert(const IndexSet& support, int start, int remaining, int index);
    virtual void dump(int level);
    virtual bool dominated(const IndexSet& b, int index1, int index2);
    virtual int index_dominated(const IndexSet& b, int index1, int index2);
    virtual int find(const IndexSet& b, int start, int remaining);
    virtual void find_diff(std::vector<int>& indices, const IndexSet& supp, int diff);
    virtual void find_diff(std::vector<int>& indices, const IndexSet& supp1, int diff1, const IndexSet& supp2, int diff2);

protected:
    int index;
};

template <class IndexSet>
class OnesTree
{
public:
    OnesTree(const std::vector<IndexSet>& supports, int num);
    OnesTree();
    ~OnesTree();

    void insert(const std::vector<IndexSet>& supports, int num);
    void insert(const IndexSet& support, int index);
    bool dominated(const IndexSet& b, int index1, int index2);
    int index_dominated(const IndexSet& b, int index1, int index2);
    int find(const IndexSet& b);
    void dump();

    void find_diff(std::vector<int>& indices, const IndexSet& supp, int diff);
    void find_diff(std::vector<int>& indices, const IndexSet& supp1, int diff1, const IndexSet& supp2, int diff2);

private:
    OnesTreeNode<IndexSet>* root;
};

} // namespace _4ti2_

#include "groebner/OnesTree.tpp"

#endif
