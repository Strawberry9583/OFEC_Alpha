#include "LinkSort.h"
//#include "../functional.h"
#include "quick_sort.h"
#include <algorithm>
#include <iostream>
#include <bitset>

namespace NDS {

	// P.S.!! data[i][j] means solution[j]'s [i]th objective !!
	void LinkSort(const std::vector<std::vector<double>>& data, std::vector<int>& rank, int& comp)
	{
		const int M = data.size();
		if (M == 0)
			return;
		const int N = data[0].size();
		std::vector<std::vector<double>> pop = data;

		std::vector<std::vector<int>> SolSeq_byOneObj(N); // each row stores indexs of solutions in ascending order by objective value
		for (int i = 0; i < N; ++i) {
			comp += quick_sort(pop, SolSeq_byOneObj[i], i);
		}
		int cur_rank = 0;
		int num_not_ranked = M; // number of solutions not ranked
		bool* nominated_last = new bool[M]; //whether in the candidate of last generation
		bool* nominated_cur = new bool[M]; // whether already in the candidate of curent generation
		std::vector<int> candidate; // the indexs of solutions nominated 
		candidate.reserve(M);
		int link;
		memset(nominated_last, true, M);
		memset(nominated_cur, false, M);
		while (num_not_ranked > 0) {
			bool flag_found = false;
			for (int i = 0; i < N; ++i) {
				for (auto j = SolSeq_byOneObj[i].begin(); j != SolSeq_byOneObj[i].end(); ++j) {
					link = *j;
					flag_found = true;
					break;
				}
				if (flag_found)
					break;
			}
			while (true) {
				for (int i = 0; i < N; ++i)
					for (auto j = SolSeq_byOneObj[i].begin(); j != SolSeq_byOneObj[i].end(); ++j) 
						if (*j == link)
							break;
						else if (nominated_last[*j] && !nominated_cur[*j]) {
							nominated_cur[*j] = true;
							candidate.push_back(*j);
						}
				rank[link] = cur_rank;
				num_not_ranked--;
				for (int i = 0; i < N; ++i)
					for (auto j = SolSeq_byOneObj[i].begin(); j != SolSeq_byOneObj[i].end(); ++j)
						if (*j == link) {
							SolSeq_byOneObj[i].erase(j);
							break;
						}
				if (candidate.size() == 0) {
					memset(nominated_last, true, M);
					break;
				}
				else if (candidate.size() == 1) {
					link = candidate.front();
					rank[link] = cur_rank;
					num_not_ranked--;
					for (int i = 0; i < N; ++i) {
						SolSeq_byOneObj[i].erase(std::remove(SolSeq_byOneObj[i].begin(), SolSeq_byOneObj[i].end(), link));
					}
					memset(nominated_last, true, M);
					nominated_cur[candidate.front()] = false;
					candidate.clear();
					break;
				}
				else {
					for (int i = 0; i < M; ++i)
						nominated_last[i] = nominated_cur[i];
					memset(nominated_cur, false, M);
					link = candidate.front();
					candidate.clear();
				}
			}
			cur_rank++;
		}
		delete[] nominated_cur;
		delete[] nominated_last;
	}
}