#include "LinkSort.h"
#include "../functional.h"

namespace NDS {

	// P.S.!! data[i][j] means solution[j]'s [i]th objective !!
	void LinkSort(std::vector<std::vector<double>>& data, std::map<int, int>& rank, int& comp)
	{
		const int N = data.size();
		const int M = data[0].size();
		std::vector<std::vector<int>> obj_seq(N); // each row stores index of solutions in sequence of each obective
		for (int i = 0; i < N; ++i)
			comp += OFEC::quick_sort(data[i], M, obj_seq[i]);
		int cur_rank = 0;
		std::vector<int> idx_ranked; // indexs of solutions ranked
		std::vector<bool> bool_ranked(M, false); // whether is ranked
		while (idx_ranked.size() < M) {
			std::vector<bool> nominated_last(M, true); //whether in the candidate of last generation
			std::list<int> candidate; // the indexs of solutions nominated 
			int link;

			bool temp_found = false;
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < M; ++j) {
					if (!bool_ranked[obj_seq[i][j]]) {
						link = obj_seq[i][j];
						temp_found = true;
						break;
					}
				}
				if (temp_found)
					break;
			}

			while (true) {
				std::vector<bool> nominated_cur(M, false); // whether already in the candidate of curent generation
				for (int i = 0; i < N; ++i) {
					for (int j = 0; j < M; ++j) {
						if (!bool_ranked[obj_seq[i][j]]) {
							if (obj_seq[i][j] == link)
								break;
							else if (nominated_last[obj_seq[i][j]]) {
								if (!nominated_cur[obj_seq[i][j]]) {
									nominated_cur[obj_seq[i][j]] = true;
									candidate.push_back(obj_seq[i][j]);
								}
							}
						}
					}
				}
				rank[link] = cur_rank;
				idx_ranked.push_back(link);
				bool_ranked[link] = true;
				if (candidate.size() == 0)
					break;
				if (candidate.size() == 1) {
					rank[candidate.front()] = cur_rank;
					idx_ranked.push_back(candidate.front());
					bool_ranked[candidate.front()] = true;
					break;
				}
				for (auto x : nominated_last)
					x = false;
				for (auto x : candidate)
					nominated_last[x] = true;
				link = candidate.front();
				candidate.clear();
			}
			cur_rank++;
		}
	}

}