#ifndef FAST_NONDOMINATED_SORT_H
#define FAST_NONDOMINATED_SORT_H

#include<map>
#include <list>

#include "../functional.h"

namespace NDS {
	class fast_sort {
	public:
		void sort() {

			for (int i = 0; i < m_popsize; i++) {
				rank[i] = -1;
			}

			auto i = cset.begin();
			for (int k = 0; k < m_popsize; k++, ++i) {
				for (int j = 0; j<m_popsize; j++) {
					if (k != j) {
						auto compare_result = OFEC::objective_compare(m_pop[j].second, m_pop[k].second, OFEC::optimization_mode::Minimization);
						m_num += compare_result.second;
						if (compare_result.first == OFEC::dominationship::Dominating) {//*m_pop[j]>*m_pop[k]
							rank_[k]++;
						}
						else if (compare_result.first == OFEC::dominationship::Dominated) {//*m_pop[k]>*m_pop[j]
							(*i)[count[k]] = j;
							count[k]++;
						}
					}
				}
			}

			int m_curRank = 0;
			std::vector<int> rank2(m_popsize);
			while (1)
			{
				int stop_count = 0;
				for (int k = 0; k<m_popsize; k++)
					rank2[k] = rank_[k];
				auto i = cset.begin();
				for (int k = 0; k<m_popsize; k++, ++i)
				{
					if (rank[k] == -1 && rank_[k] == 0)
					{
						rank[k] = m_curRank;
						for (int j = 0; j<count[k]; j++)
						{
							int id = (*i)[j];
							rank2[id]--;
							stop_count++;
						}
					}
				}

				for (int k = 0; k<m_popsize; k++)
					rank_[k] = rank2[k];

				m_curRank++;
				if (stop_count == 0)
					break;
			}
		}
		std::map<int, int>& ranking() {

			for (int i = 0; i < m_popsize; ++i) {
				m_label.emplace(m_pop[i].first, rank[i]);
			}
			return m_label;
		}
		int number() {
			return m_num;
		}
		fast_sort(std::vector<std::pair<int, std::vector<double>>> && data) :m_popsize(data.size()), m_pop(std::move(data)), rank_(m_popsize), rank(m_popsize), count(m_popsize),
			cset(m_popsize, std::vector<int>(m_popsize)) {
			if (!data.empty())
				m_objsnum = data[0].second.size();
		}
		fast_sort(const std::vector<std::pair<int, std::vector<double>>> & data) :m_popsize(data.size()), m_pop(data), rank_(m_popsize), rank(m_popsize), count(m_popsize),
			cset(m_popsize, std::vector<int>(m_popsize)) {
			if (!data.empty())
				m_objsnum = data[0].second.size();
		}
	private:
		int m_objsnum;
		int m_popsize;
		std::vector<std::pair<int, std::vector<double>>> m_pop;
		std::vector<int> rank_, count, rank;
		std::list<std::vector<int> > cset;
		int m_num = 0;

		std::map<int, int> m_label;
	};
}




#endif // !FAST_NONDOMINATED_SORT_H