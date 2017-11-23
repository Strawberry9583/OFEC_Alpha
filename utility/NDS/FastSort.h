#ifndef FAST_NONDOMINATED_SORT_H
#define FAST_NONDOMINATED_SORT_H

#include<map>
#include <list>

#include "../functional.h"

namespace NDS {
	class fast_sort {
	public:
		fast_sort(std::vector<std::pair<int, std::vector<double>>> && data);
		fast_sort(const std::vector<std::pair<int, std::vector<double>>> & data);
		void sort();
		std::vector<int>& rank_result() { return rank; }
		int number() { return m_num; }
		int rank_num() { return m_rank_num; }
	private:
		size_t m_objsnum;
		size_t m_popsize;
		std::vector<std::pair<int, std::vector<double>>> m_pop;
		std::vector<int> rank_, count, rank;
		std::list<std::vector<int> > cset;
		int m_num = 0;
		int m_rank_num;

		std::map<int, int> m_label;
	};
}




#endif // !FAST_NONDOMINATED_SORT_H