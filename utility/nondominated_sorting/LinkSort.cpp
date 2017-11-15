#include "LinkSort.h"
#include "../functional.h"
#include <iostream>
#include <time.h>

namespace NDS {

	// P.S.!! data[i][j] means solution[j]'s [i]th objective !!
	void LinkSort(std::vector<std::vector<double>>& data, std::vector<int>& rank, int& comp)
	{
		clock_t time_start, time_finish, time_part1(0), time_part2(0), time_part3(0), time_part4(0);
		const int N = data.size();
		const int M = data[0].size();
		std::vector<std::list<int>*> obj_seq(N); 
		//std::vector<std::vector<int>*> obj_seq(N);// each row stores indexs of solutions in sequence of each obective
		for (int i = 0; i < N; ++i) {
			std::vector<int> temp_row_obj_seq;
			//part1 start
			time_start = clock();
			comp += OFEC::quick_sort(data[i], M, temp_row_obj_seq);
			obj_seq[i] = new std::list<int>(temp_row_obj_seq.begin(), temp_row_obj_seq.end());
			time_finish = clock();
			time_part1 += (time_finish - time_start);
			//part1 finish
			//obj_seq[i] = new std::vector<int>;
			//obj_seq[i]->reserve(M);
			//comp += OFEC::quick_sort(data[i], M, *obj_seq[i]);
		}
		int cur_rank = 0;
		int num_not_ranked = M; // number of solutions not ranked
		std::vector<bool> nominated_last(M, true); //whether in the candidate of last generation
		std::vector<bool> nominated_cur(M, false); // whether already in the candidate of curent generation
		std::vector<int> candidate; // the indexs of solutions nominated 
		candidate.reserve(M);
		int link;
		while (num_not_ranked > 0) {
			//part2 start
			time_start = clock();
			bool flag_found = false;
			for (int i = 0; i < N; ++i) {
				for (auto j = obj_seq[i]->begin(); j != obj_seq[i]->end(); ++j) {
					link = *j;
					flag_found = true;
					break;
				}
				if (flag_found)
					break;
			}
			time_finish = clock();
			time_part2 += (time_finish - time_start);
			//part2 finish
			while (true) {
				//part3 start
				time_start = clock();
				for (int i = 0; i < N; ++i) {
					for (auto j = obj_seq[i]->begin(); j != obj_seq[i]->end(); ++j) {
						if (*j == link)
							break;
						if (nominated_last[*j]) {
							if (!nominated_cur[*j]) {
								nominated_cur[*j] = true;
								candidate.push_back(*j);
							}
						}
					}
				}
				time_finish = clock();
				time_part3 += (time_finish - time_start);
				//part3 finish
				//part4 start
				time_start = clock();
				rank[link] = cur_rank;
				num_not_ranked--;
				for (int i = 0; i < N; ++i) {
					for (auto j = obj_seq[i]->begin(); j != obj_seq[i]->end(); ++j)
						if (*j == link) {
							obj_seq[i]->erase(j);
							break;
						}
				}
				time_finish = clock();
				time_part4 += (time_finish - time_start);
				//part4 finish
				if (candidate.size() == 0) {
					for (auto& x : nominated_last)
						x = true;
					candidate.clear();
					break;
				}
				else if (candidate.size() == 1) {
					link = candidate.front();
					rank[link] = cur_rank;
					num_not_ranked--;
					for (int i = 0; i < N; ++i) {
						for (auto j = obj_seq[i]->begin(); j != obj_seq[i]->end(); ++j)
							if (*j == link) {
								obj_seq[i]->erase(j);
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
		for (int i = 0; i < N; ++i)
			delete obj_seq[i];
		std::cout << "time_part1 cost: " << time_part1 << " milliseconds." << std::endl;
		std::cout << "time_part2 cost: " << time_part2 << " milliseconds." << std::endl;
		std::cout << "time_part3 cost: " << time_part3 << " milliseconds." << std::endl;
		std::cout << "time_part4 cost: " << time_part4 << " milliseconds." << std::endl;
	}

}