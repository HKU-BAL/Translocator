/*
 * IntervallTree.h
 *
 *  Created on: Jun 23, 2015
 *      Author: fsedlaze
 */

#ifndef TREE_INTERVALLTREE_H_
#define TREE_INTERVALLTREE_H_

#include <vector>

#include "TNode.h"
#include "IntervallContainer.h"



class IntervallTree:public IntervallContainer {
private:
	int max(int, int);
	TNode * srl(TNode *&);
	TNode * drl(TNode *&);
	TNode * srr(TNode *&);
	TNode * drr(TNode *&);
	void careful_screening(Breakpoint *& new_break, TNode *p);
public:
	void insert(Breakpoint * point, TNode *&);
	void insert_ref(Breakpoint * point, TNode *&);

	void insert_existant(Breakpoint * new_break, TNode *&p);
	void del(Breakpoint * point, TNode *&);
	int deletemin(TNode *&);
	void find(Breakpoint * point, TNode *&);
	bool overlaps(long start, long stop,TNode *p);
	TNode * findmin(TNode*);
	TNode * findmax(TNode*);
	void clear(TNode *&);
	void copy(TNode * &, TNode *&);
	TNode * nodecopy(TNode *&);
	void preorder(TNode*);
	void inorder(TNode * p);
	void postorder(TNode*);
	int bsheight(TNode*);
	void get_breakpoints(TNode *p,std::vector<Breakpoint *> & points);
	int nonodes(TNode*);
	void collapse_intervalls(TNode *&p);
	void print(TNode *p);
};

#endif /* TREE_INTERVALLTREE_H_ */
