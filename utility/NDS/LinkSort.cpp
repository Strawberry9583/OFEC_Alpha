#include "LinkSort.h"
#include "../functional.h"
#include "quick_sort.h"
#include <algorithm>
#include <iostream>
#include <limits.h>
#include <string.h>
#ifdef USING_CONCURRENT
#include <thread>
#endif // USING_CONCURRENT
#include <chrono>
namespace NDS {
	void LinkSort(const std::vector<std::vector<double>>& data, std::vector<int>& rank, int& comp){
		std::chrono::time_point<std::chrono::system_clock> start_time;
		std::chrono::milliseconds time_cost;
		time_cost = time_cost.zero();
		const int M = data.size(); // number of solution
		if (M == 0)
			return;
		const int N = data[0].size(); // number of obective
		if (rank.empty() || rank.size() != M)
			rank.resize(M);

		// ******************************************

		std::vector<std::vector<int>> SeqByObj(N);
        #ifdef USING_CONCURRENT
		int TaskSize = N;
		int numTask = std::thread::hardware_concurrency();
		if (numTask > TaskSize) numTask = TaskSize;
		std::vector<std::thread> thrd;
		for (int i = 0; i < numTask; ++i) {
			std::vector<int> ObjIdxs;
			for (int idx = i; idx < TaskSize; idx += numTask)
				ObjIdxs.push_back(idx);
			thrd.push_back(std::thread(ParallelQuickSort, std::cref(data), std::ref(SeqByObj), std::move(ObjIdxs)));
		}
		for (auto&t : thrd) 
			t.join();
        #else
		for (int i = 0; i < N; ++i)
			comp += quick_sort(data, SeqByObj[i], i);
        #endif // USING_CONCURRENT
		std::vector<LS_list> SeqByObj_Lists(N); // same as SeqByObj but in form of LS_list
		std::vector<std::vector<LS_node*>> PosInObjLists(M); // PosInObjLists[i] stores solution[i]'s all LS_node addresses 
		for (int i = 0; i < N; ++i) {
			for (int SolIdx : SeqByObj[i])
				PosInObjLists[SolIdx].push_back(SeqByObj_Lists[i].push_back(SolIdx));
		}

		std::vector<std::vector<int>> SolStas(M); // SolStas[i] stores solution[i]'s all single objective sequence numbers
		for (int i = 0; i < N; ++i) {
			int SeqNum = 0;
			for (LS_node* iter = SeqByObj_Lists[i].begin(); iter != nullptr; iter = iter->m_next) {
				SolStas[iter->m_value].push_back(SeqNum);
				SeqNum++;
			}
		}

		//std::vector<int> MinIdxs(M); // MinIdxs[i] means index of solution[i]'s minimum single objective sequence number
		//std::vector<int> SumVals(M); // SumVals[i] means sum of solution[i]'s all single objective sequence numbers
		//for (int i = 0; i < M; ++i) {
		//	int min_val = INT_MAX; // value of solution's minimum single objective sequence number
		//	int min_idx; // index of solution's minimum single objective sequence number
		//	for (int ObjIdx = 0; ObjIdx < N; ++ObjIdx) {
		//		if (SolStas[i][ObjIdx] < min_val) {
		//			min_val = SolStas[i][ObjIdx];
		//			min_idx = ObjIdx;
		//		}
		//		SumVals[i] += SolStas[i][ObjIdx];
		//	}
		//	MinIdxs[i] = min_idx;
		//}

		std::vector<int> MaxIdxs(M);
		std::vector<int> MinIdxs(M); // MinIdxs[i] means index of solution[i]'s minimum single objective sequence number
		std::vector<int> SumVals(M); // SumVals[i] means sum of solution[i]'s all single objective sequence numbers
		for (int i = 0; i < M; ++i) {
			int max_val = -1;
			int max_idx;
			int min_val = INT_MAX; // value of solution's minimum single objective sequence number
			int min_idx; // index of solution's minimum single objective sequence number
			for (int ObjIdx = 0; ObjIdx < N; ++ObjIdx) {
				if (SolStas[i][ObjIdx] < min_val) {
					min_val = SolStas[i][ObjIdx];
					min_idx = ObjIdx;
				}
				if (SolStas[i][ObjIdx] > max_val) {
					max_val = SolStas[i][ObjIdx];
					max_idx = ObjIdx;
				}
				SumVals[i] += SolStas[i][ObjIdx];
			}
			MinIdxs[i] = min_idx;
			MaxIdxs[i] = max_idx;
		}

		//std::vector<int> SumVals(M); // SumVals[i] means sum of solution[i]'s all single objective sequence numbers
		//for (int i = 0; i < M; ++i)
		//	for (int ObjIdx = 0; ObjIdx < N; ++ObjIdx) 
		//		SumVals[i] += SolStas[i][ObjIdx];
		//std::vector<std::vector<int>> MaxIdxss(M); //MinIdxs[i][j] means index of solution[i]'s [j]th maximum single objective sequence number
		//for (int i = 0; i < M; ++i)
		//	OFEC::quick_sort(SolStas[i], N, MaxIdxss[i]);

		std::vector<int> SeqBySumVals(M); // sequence of solution sorted by sum of all single objective sequence numbers
		OFEC::quick_sort(SumVals, M, SeqBySumVals);
		LS_list SeqBySumVals_Lists; // Same as SeqByMinVals but in form of LS_list
		for (const int SolIndex : SeqBySumVals)
			PosInObjLists[SolIndex].push_back(SeqBySumVals_Lists.push_back(SolIndex));

		// ******************************************

		int CurRankNumber = 0; // current rank number
		int num_not_ranked = M; // number of solutions not ranked
		std::vector<int> CurRankCandidate; // the candidate solutions of current rank
		CurRankCandidate.reserve(M);
		bool* InCurRankCandiate = new bool[M]; // whether in CurRankCandidate
		memset(InCurRankCandiate, false, M);
		int link; // the solution chosen to generate CurRankCandidate
		while (num_not_ranked > 0) {
			// set link
			link = SeqBySumVals_Lists.begin()->m_value;
			rank[link] = CurRankNumber;
			num_not_ranked--;
			// generate CurRankCandidate by link
			for (auto SeqByObj_List : SeqByObj_Lists) {
				for (auto iter = SeqByObj_List.begin(); iter != nullptr; iter = iter->m_next) {
					if (iter->m_value == link)
						break;
					else if (!InCurRankCandiate[iter->m_value]) {
						InCurRankCandiate[iter->m_value] = true;
						CurRankCandidate.push_back(iter->m_value);
					}
				}
			}
			// remove solution[link]'s all LS_nodes
			for (int i = 0; i < N; ++i)
				SeqByObj_Lists[i].erase(PosInObjLists[link][i]);
			SeqBySumVals_Lists.erase(PosInObjLists[link][N]);
			// filter the CurRankCandidate
#ifdef USING_CONCURRENT
			numTask = std::thread::hardware_concurrency();
			TaskSize = CurRankCandidate.size();
			if (numTask > TaskSize) numTask = TaskSize;
			thrd.clear();
			for (int i = 0; i < numTask; ++i) {
				std::vector<int> candidates;
				for (int idx = i; idx < TaskSize; idx += numTask)
					candidates.push_back(CurRankCandidate[idx]);
				thrd.push_back(std::thread(ParallelFilter, std::move(candidates), std::ref(SeqByObj_Lists), std::cref(MinIdxs), N, std::cref(SolStas), InCurRankCandiate));
			}
			for (auto&t : thrd) 
				t.join();
#else
			start_time = std::chrono::system_clock::now();
			for (auto candidate : CurRankCandidate) {
				bool FlagInCurRank(true); // whether candidate is in current rank 
				for (auto iter = SeqByObj_Lists[MinIdxs[candidate]].begin(); iter != nullptr; iter = iter->m_next) {
				//for (auto iter = SeqByObj_Lists[MaxIdxss[candidate][N - 1]].begin(); iter != nullptr; iter = iter->m_next) {
					if (iter->m_value == candidate)
						break;
					else {
						if (SolStas[iter->m_value][MaxIdxs[iter->m_value]] > SolStas[candidate][MaxIdxs[iter->m_value]])
							continue;
						// check whether solution[iter->m_value] donminate solution[candidate] or not
						bool FlagDominate(true);
						for (int i = 0; i < N; ++i) {
							if (i != MaxIdxs[iter->m_value] && SolStas[iter->m_value][i] > SolStas[candidate][i]) {
							//if (SolStas[iter->m_value][i] > SolStas[candidate][i]) {
							//int objidx = MaxIdxss[iter->m_value][i];
							//if (SolStas[iter->m_value][objidx] > SolStas[candidate][objidx]) {
								FlagDominate = false;
								break;
							}
						}
						if (FlagDominate) {
							FlagInCurRank = false;
							break;
						}
					}
				}
				if (!FlagInCurRank)
					InCurRankCandiate[candidate] = false;
			}
			time_cost += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
#endif // USING_CONCURRENT
			// set rank for filtered CurRankCandidate and remove their LS_nodes
			for (auto candidate : CurRankCandidate) {
				if (InCurRankCandiate[candidate]) {
					rank[candidate] = CurRankNumber;
					num_not_ranked--;
					for (int i = 0; i < N; ++i)
						SeqByObj_Lists[i].erase(PosInObjLists[candidate][i]);
					SeqBySumVals_Lists.erase(PosInObjLists[candidate][N]);
				}
			}
			// goto next rank
			CurRankNumber++;
			CurRankCandidate.clear();
			memset(InCurRankCandiate, false, M);
		}
		delete InCurRankCandiate;
		std::cout << "Dominance check cost:" << time_cost.count() << std::endl;
	}
#ifdef USING_CONCURRENT
	void ParallelFilter(const std::vector<int>&& candidates, std::vector<LS_list>& SeqByObj_Lists, const std::vector<int>& MinIdxs, const int N, const std::vector<std::vector<int>>& SolStas, bool* InCurRankCandiate) {
		for (int candidate : candidates) {
			bool FlagInCurRank(true); // whether candidate is in current rank 
			for (auto iter = SeqByObj_Lists[MinIdxs[candidate]].begin(); iter != nullptr; iter = iter->m_next) {
				if (iter->m_value == candidate)
					break;
				else {
					// check whether solution[iter->m_value] donminate solution[candidate] or not
					bool FlagDominate(true);
					for (int i = 0; i < N; ++i)
						if (SolStas[iter->m_value][i] > SolStas[candidate][i]) {
							FlagDominate = false;
							break;
						}
					if (FlagDominate) {
						FlagInCurRank = false;
						break;
					}
				}
			}
			if (!FlagInCurRank)
				InCurRankCandiate[candidate] = false;
		}
	}
	void ParallelQuickSort(const std::vector<std::vector<double>>& data, std::vector<std::vector<int>>& SeqByObj, const std::vector<int>&& ObjIdxs) {
		for (int ObjIdx : ObjIdxs)
			quick_sort(data, SeqByObj[ObjIdx], ObjIdx);
	}
#endif // USING_CONCURRENT
}