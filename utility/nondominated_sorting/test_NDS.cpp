#include "../../utility/nondominated_sorting/FastNondominatedSort.h"
#include "../../utility/nondominated_sorting/CornerSort.h"
#include "../../utility/nondominated_sorting/T_ENS.h"
#include "../../utility/nondominated_sorting/LinkSort.h"
#include "../../utility/nondominated_sorting/preSorting.h"
#include "../../utility/nondominated_sorting/nodes_initialize.h"

#include <time.h>
#include <fstream>

int main() {
	const int data_size(10000);
	const int obj_num(3);
	OFEC::random rand1(0.2);
	clock_t start(0), end(0);

	NDS::circle_distribution d1(obj_num, data_size, rand1);
	std::vector<std::vector<double>> data = d1.get_data();

	/*for (int i = 0; i < 500; ++i)
		data[i][0] = 0;*/

	std::vector<std::vector<double>> data_T(obj_num, std::vector<double>(data_size));
	for (int i = 0; i < data_size; ++i)
		for (int j = 0; j < obj_num; ++j)
			data_T[j][i] = data[i][j];

	std::vector<int> no;
	std::vector<std::pair<int, std::vector<double>>> md;
	for (int i = 0; i < data_size; ++i) {
		no.push_back(i);
		md.emplace_back(no[i], data[i]);
	}

	const int run_num = 10;
	clock_t time_cost(0);

	std::vector<int> ls_rank(data_size);
	std::vector<int> te_rank(data_size);
	std::vector<int> cs_rank(data_size);

	for (int runID = 0; runID < run_num; ++runID) {
		int ls_com = 0;
		start = clock();
		NDS::LinkSort(data_T, ls_rank, ls_com);
		end = clock();
		std::cout << end - start << " ";
		time_cost += end - start;
	}
	std::cout << "\nLinkSort cost:" << time_cost / run_num << " milliseconds\n\n";
	time_cost = 0;

	for (int runID = 0; runID < run_num; ++runID) {
		int te_com = 0;
		start = clock();
		NDS::T_ENS(data, te_com, te_rank);
		end = clock();
		std::cout << end - start << " ";
		time_cost += end - start;
	}
	std::cout << "\nT-ENS cost:" << time_cost / run_num << " milliseconds\n\n";
	time_cost = 0;

	for (int runID = 0; runID < run_num; ++runID) {
		double** POP = new double*[data_size];
		for (int i = 0; i < data_size; ++i)
			POP[i] = data[i].data();
		std::vector<int> cs_comp(data_size, 0);
		int num_comp;
		start = clock();
		NDS::cornerSort(POP, obj_num, data_size, cs_rank.data(), cs_comp.data(), num_comp);
		end = clock();
		std::cout << end - start << " ";
		time_cost += end - start;
	}
	std::cout << "\nCorner_Sort cost:" << time_cost / run_num << " milliseconds\n\n";
	time_cost = 0;

	NDS::fast_sort fs(md);
	start = clock();
	fs.sort();
	end = clock();
	std::cout << "\nFNS cost:" << end - start << " milliseconds\n\n";
	std::vector<int>& fs_rank = fs.rank_result();

	std::cout << (fs_rank == te_rank ? "fs_rank == te_rank" : "fs_rank != te_rank") << std::endl;
	std::cout << (fs_rank == ls_rank ? "fs_rank == ls_rank" : "fs_rank != ls_rank") << std::endl;
	std::cout << (fs_rank == cs_rank ? "fs_rank == cs_rank" : "fs_rank != cs_rank") << std::endl;

	system("pause");
	return 0;
}