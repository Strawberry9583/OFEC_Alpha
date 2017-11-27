#include "DeductiveSort.h"
#include <string.h>
#include "../functional.h"

namespace NDS {
	void DeductiveSort(const std::vector<std::vector<double>>& data, std::vector<int>& rank, int & comp){
		const int N = data.size(); // size of data
		if (N == 0) return;
		const int M = data.front().size(); // number of objective
		int x = 0; // front number
		int f = 0; // number sorted
		bool* F = new bool[N]; // ranked flag
		memset(F, false, N);
		bool* D = new bool[N]; // dominated flag
		while (f < N) {
			memset(D, false, N);
			for (int i = 0; i < N; ++i) {
				if (!D[i] && !F[i]) {
					for (int j = i + 1; j < N; ++j) {
						if (!D[j] && !F[j]) {
							std::pair<OFEC::dominationship, int> result = OFEC::objective_compare(data[j], data[i], OFEC::optimization_mode::Minimization);
							OFEC::dominationship d = result.first;
							comp += result.second;
							if (d == OFEC::dominationship::Dominated)
								D[j] = true;
							else if (d == OFEC::dominationship::Dominating) {
								D[i] = true;
								break;
							}
						}
					}
					if (!D[i]) {
						rank[i] = x;
						f++;
						F[i] = true;
					}
				}
			}
			x++;
		}
		delete D;
		delete F;
	}
}