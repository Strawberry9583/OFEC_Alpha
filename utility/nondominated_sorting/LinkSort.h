#ifndef LINKSORT_H
#define LINKSORT_H

#include <vector>
#include <map>
#include <list>

namespace NDS {
	void LinkSort(std::vector<std::vector<double>>& data, std::map<int, int>& rank, int& comp);
}

#endif // !LINKSORT_H
