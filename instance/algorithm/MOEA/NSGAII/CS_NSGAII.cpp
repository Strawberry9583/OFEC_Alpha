#include "CS_NSGAII.h"

namespace OFEC {
	void CS_NSGAII::sort() {
		const size_t data_size(m_offspring.size());
		int obj_num = global::ms_global->m_problem->objective_size();
		double** POP = new double*[data_size];
		for (size_t i = 0; i < data_size; ++i)
			POP[i] = m_offspring[i].get_objective().data();
		int* cs_rank = new int[data_size];
		int* cs_com = new int[obj_num] { 0 };
		NDS::cornerSort(POP, obj_num, data_size, cs_rank, cs_com, m_objcomp);
		for (size_t i = 0; i < data_size; ++i)
			m_offspring[i].set_rank(cs_rank[i]);
		delete[] POP;
		delete[] cs_rank;
		delete[] cs_com;
	}
}
