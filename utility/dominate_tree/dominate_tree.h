#ifndef DOMINATE_TREE_H
#define DOMINATE_TREE_H

#include <vector>
#include "../core/definition.h"

using namespace std;

namespace dominate_tree {
	class node {
	private:
		size_t ranking;
		vector<double> m_obj;
		vector<node*> children;
	public:
		void add_solution() {

		}
		void delete_solution() {

		}


	};
}
#endif // !DOMINATE_TREE_H
