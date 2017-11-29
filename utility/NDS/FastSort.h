#ifndef FAST_NONDOMINATED_SORT_H
#define FAST_NONDOMINATED_SORT_H

#include<map>
#include <list>
#include <vector>
#include "../functional.h"

namespace NDS {
	int FastSort(const std::vector<std::vector<double>>& data, std::vector<int>& rank, int& num_comp);
#ifdef USING_CONCURRENT
	void ParallelCompare(int popsize, const std::vector<int>&& ks, const std::vector<std::vector<double>>& data, std::vector<int>& rank_, std::vector<int>& count, std::vector<std::vector<int>>& cset);
#endif // USING_CONCURRENT
}




#endif // !FAST_NONDOMINATED_SORT_H