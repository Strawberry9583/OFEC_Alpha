#ifndef DEDUCTIVESORT_H
#define DEDUCTIVESORT_H

/*K.M.Clymont and E.Keedwell, ��Deductive sort and climbing sort :
new methods for non - dominated sorting, �� Evolutionary Computation,
vol. 20, no. 1, pp. 1�C26, 2012.*/

#include <vector>

namespace NDS {
	void DeductiveSort(const std::vector<std::vector<double>>& data, std::vector<int>& rank, int& comp);
}

#endif // !DEDUCTIVESORT_H