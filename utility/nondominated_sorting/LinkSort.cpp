#include "LinkSort.h"
#include "../functional.h"
#include <iostream>

namespace NDS {

	// P.S.!! data[i][j] means solution[j]'s [i]th objective !!
	void LinkSort(const std::vector<std::vector<double>>& data, std::vector<int>& rank, int& comp)
	{
		const int N = data.size();
		const int M = data[0].size();
		//std::vector<std::list<int>> SolSeq_byOneObj(N); // each row stores indexs of solutions in ascending order by objective value
		std::vector<std::vector<int>> SolSeq_byOneObj(N); // each row stores indexs of solutions in ascending order by objective value
		//std::vector<int> ObjSeqSum_byOneObj(M); // each row stores sequences of solutions' one objectives in ascending order by solution index, last row stores sum of each column
		//std::vector<int> SolSeq_byAllObj(M);  //stores indexs of solutions in ascending order by all objective value
		for (int i = 0; i < N; ++i) {
			//std::vector<int> temp_row;
			//comp += OFEC::quick_sort(data[i], M, temp_row);
			//SolSeq_byOneObj[i] = std::list<int>(temp_row.begin(), temp_row.end());
			comp += OFEC::quick_sort(data[i], M, SolSeq_byOneObj[i]);
			//for (int j = 0; j < M; ++j)
			//	ObjSeqSum_byOneObj[SolSeq_byOneObj[i][j]] += j;
		}
		//OFEC::quick_sort(ObjSeqSum_byOneObj, M, SolSeq_byAllObj);
		int cur_rank = 0;
		int num_not_ranked = M; // number of solutions not ranked
		std::vector<bool> nominated_last(M, true); //whether in the candidate of last generation
		std::vector<bool> nominated_cur(M, false); // whether already in the candidate of curent generation
		std::vector<int> candidate; // the indexs of solutions nominated 
		candidate.reserve(M);
		int link;
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
			//link = SolSeq_byAllObj.front();
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
		/*		for (auto j = SolSeq_byAllObj.begin(); j != SolSeq_byAllObj.end(); ++j)
					if (*j == link) {
						SolSeq_byAllObj.erase(j);
						break;
					}*/
				if (candidate.size() == 0) {
					for (auto& x : nominated_last)
						x = true;
					break;
				}
				else if (candidate.size() == 1) {
					link = candidate.front();
					rank[link] = cur_rank;
					num_not_ranked--;
					for (int i = 0; i < N; ++i) {
						for (auto j = SolSeq_byOneObj[i].begin(); j != SolSeq_byOneObj[i].end(); ++j)
							if (*j == link) {
								SolSeq_byOneObj[i].erase(j);
								break;
							}
					}
					for (auto& x : nominated_last)
						x = true;
					nominated_cur[candidate.front()] = false;
					candidate.clear();
					break;
				}
				else {
					nominated_last = nominated_cur;
					for (auto x : candidate)
						nominated_cur[x] = false;
					link = candidate.front();
					candidate.clear();
				}
			}
			cur_rank++;
		}
	}

}