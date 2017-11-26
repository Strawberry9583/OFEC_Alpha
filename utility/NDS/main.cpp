#include "FastSort.h"
#include "CornerSort.h"
#include "T_ENS.h"
#include "LinkSort.h"
#include "static_population.h"

#include <time.h>
#include <fstream>

int main(int argc, char* argv[]) {
	int data_size(atoi(argv[1]));
	int obj_num(atoi(argv[2]));
	int rank_num(atoi(argv[3]));
	const int run_num(10);
	std::cout << "data_size:" << data_size << "  obj_num:" << obj_num << "  rank_num:" << rank_num << "  run_num:" << run_num << std::endl;
	clock_t start(0), end(0), time_cost(0);
	NDS::uniform_population u1(obj_num, data_size, 0.5);
	std::vector<std::vector<double>> data = u1.generate_new(rank_num);
	std::cout << "done" << std::endl;

	//for (int i = 0; i < 500; ++i)
	//	data[i][1] = 0;

	//for (int i = 800; i < 1800; ++i)
	//	data[i][1] = 1;

	//for (int i = 2000; i < 2800; ++i)
	//	data[i][0] = 1;

	std::vector<int> no;
	std::vector<std::pair<int, std::vector<double>>> md;
	for (int i = 0; i < data_size; ++i) {
		no.push_back(i);
		md.emplace_back(no[i], data[i]);
	}

	std::vector<int> ls_rank(data_size);
	std::vector<int> te_rank(data_size);
	std::vector<int> cs_rank(data_size);
	std::vector<int> lg_rank(data_size);
	std::vector<int> fs_rank(data_size);

	for (int runID = 0; runID < run_num; ++runID) {
		int ls_com = 0;
		start = clock();
		NDS::LinkSort(data, ls_rank, ls_com);
		end = clock();
		time_cost += end - start;
	}
	std::cout << "\nLS:\t" << time_cost / run_num << " milliseconds\n";
	time_cost = 0;

	for (int runID = 0; runID < run_num; ++runID) {
		int te_com = 0;
		start = clock();
		NDS::T_ENS(data, te_com, te_rank);
		end = clock();
		time_cost += end - start;
	}
	std::cout << "\nT-ENS:\t" << time_cost / run_num << " milliseconds\n";
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
		time_cost += end - start;
	}
	std::cout << "\nCS:\t" << time_cost / run_num << " milliseconds\n";
	time_cost = 0;

	//for (int runID = 0; runID < run_num; ++runID) {
	//	start = clock();
	//	NDS::ObjecitveLink_Graph_Sort lg(data_size, data);
	//	lg.initGraphLinkSort();
	//	end = clock();
	//	time_cost += end - start;
	//	lg.getRank(lg_rank);
	//}
	//std::cout << "\nLGS:\t" << time_cost / run_num << " milliseconds\n";
	//time_cost = 0;

	for (int runID = 0; runID < run_num; ++runID) {
		NDS::fast_sort fs(md);
		start = clock();
		fs.sort();
		end = clock();
		time_cost += end - start;
		fs_rank = fs.rank_result();
		rank_num = fs.rank_num();
	}
	std::cout << "\nFS:\t" << time_cost / run_num << " milliseconds\n" << std::endl;
	time_cost = 0;

	std::cout << "rank_num:\t" << rank_num << std::endl << std::endl;

	std::cout << (fs_rank == te_rank ? "fs_rank == te_rank" : "fs_rank != te_rank") << std::endl;
	std::cout << (fs_rank == ls_rank ? "fs_rank == ls_rank" : "fs_rank != ls_rank") << std::endl;
	std::cout << (fs_rank == cs_rank ? "fs_rank == cs_rank" : "fs_rank != cs_rank") << std::endl;
	//std::cout << (fs_rank == lg_rank ? "fs_rank == lg_rank" : "fs_rank != lg_rank") << std::endl;

	system("pause");
	return 0;
}