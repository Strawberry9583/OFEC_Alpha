#include "new_static_population.h"
#include <list>
#include <fstream>
#include <cmath>

namespace NDS {
	std::vector<std::vector<double>> new_uniform_population::generate_new(const int rank_num)
	{
		OFEC::uniform rand(0.5);
		std::vector<std::vector<double>> data;
		data.reserve(m_pop_size);
		int cur_rank = 1;
		int num_node = 0;
		std::vector<std::vector<double>> last_bounds;
		std::vector<std::vector<double>> cur_bounds;
		last_bounds.push_back(std::vector<double>(m_obj_num, 0.));
		while (num_node < m_pop_size) {
			std::vector<double> new_node(m_obj_num);
			double radius = cur_rank * 1. / rank_num;
			std::vector<double>& rand_bound = last_bounds[rand.next_non_standard(static_cast<size_t>(0), last_bounds.size())];
			std::vector<double> upper_bound(m_obj_num);
			for (int i = 0; i < m_obj_num; ++i) {
				double pow_sum(0);
				for (int j = 0; j < m_obj_num; ++j) {
					if (j != i)
						pow_sum += std::pow(rand_bound[j],2);
				}
				upper_bound[i] = std::sqrt(1 - pow_sum);
			}
			for (int i = 0; i < m_obj_num; ++i)
				new_node[i] = rand.next_non_standard(rand_bound[i], upper_bound[i]);
			double lenth(0);
			for (double x : new_node)
				lenth += std::pow(x, 2);
			lenth = std::sqrt(lenth);
			for (double& x : new_node)
				x = x / lenth * radius;
			data.push_back(new_node);
			cur_bounds.push_back(new_node);
			num_node++;
			if (cur_rank < rank_num && num_node == cur_rank * m_pop_size / rank_num) {
				cur_rank++;
				last_bounds = cur_bounds;
				cur_bounds.clear();
			}
		}
		return data;
	}
	std::vector<std::vector<double>> new_uniform_population::generate_output(const int rank_num) {
		std::vector<std::vector<double>> data = generate_new(rank_num);
		std::ofstream outputfile("E://new_uniform_population.csv");
		for (auto& row : data) {
			for (int i = 0; i < m_obj_num - 1; ++i)
				outputfile << row[i] << ", ";
			outputfile << row[m_obj_num - 1] << std::endl;
		}
		outputfile.close();
		return data;
	}
}