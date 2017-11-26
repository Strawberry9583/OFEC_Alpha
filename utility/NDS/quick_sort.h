#ifndef NDS_QUICK_SORT_H
#define NDS_QUICK_SORT_H

#include "../../core/algorithm/encoding.h"
#include "../../core/definition.h"
#include <algorithm>
//#include <numeric>

namespace NDS {

	template<typename T>
	bool greater(const std::vector<T>& a, const std::vector<T>& b, const int obj_idx, int& Noc) {
		size_t obj_num = a.size();
		for (size_t i = 0; i < obj_num; ++i) {
			int idx = (i + obj_idx) % obj_num;
			Noc++;
			if (a[idx] > b[idx])
				return true;
			else if (a[idx] < b[idx])
				return false;
		}
		return true;
	}

	template<typename T>
	int Partition(const std::vector<std::vector<T>>& data, std::vector<int>& A, const int obj_idx, int low, int high, int& Noc)
	{
		int pivot = A[low];
		while (low < high)
		{
			while (low < high && greater(data[A[high]], data[pivot], obj_idx, Noc))
				--high;

			A[low] = A[high]; //将比枢轴值小的元素移到左边  

			while (low < high && greater(data[pivot], data[A[low]], obj_idx, Noc))
				++low;

			A[high] = A[low]; //将比枢轴值小的元素移到右边  
		}
		A[low] = pivot;  //将枢轴值元素置于最终位置  
		return low;
	}

	template<typename T>
	void QuickSort(const std::vector<std::vector<T>>& data, std::vector<int>& A, const int obj_idx, int low, int high, int& Noc)
	{
		if (low < high) //递归跳出条件
		{
			//Partition就是上面的划分操作
			int pivot = Partition(data, A, obj_idx, low, high, Noc); //划分

			QuickSort(data, A, obj_idx, low, pivot - 1, Noc); //左半部分递归

			QuickSort(data, A, obj_idx, pivot + 1, high, Noc); //右半部分递归
		}
	}

	template<typename T>
	int quick_sort(const std::vector<std::vector<T>>& data, std::vector<int>& index, const int obj_idx = 0, bool ascending = true) {
		const int N = data.size();
		if (index.size() == 0 || index.size() != N)		index.resize(N);
		for (auto i = index.begin(); i != index.end(); ++i) *i = i - index.begin();
		int Noc(0);
		QuickSort(data, index, obj_idx, 0, index.size() - 1, Noc);
		if (!ascending) {
			for (size_t i = 0; i < data.size() / 2; i++)
				std::swap(index[i], index[N - 1 - i]);
		}
		return Noc;
	}
}

#endif // !NDS_QUICK_SORT_H

