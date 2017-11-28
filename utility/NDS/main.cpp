#include "FastSort.h"
#include "CornerSort.h"
#include "T_ENS.h"
#include "LinkSort.h"
#include "YiyaSort.h"
#include "DeductiveSort.h"
#include "static_population.h"
#include "new_static_population.h"
#include "nodes_initialize.h"

#include <chrono>
//#include <time.h>
#include <fstream>

int main(int argc, char* argv[]) {
	int data_size(atoi(argv[1]));
	int obj_num(atoi(argv[2]));
	//int rank_num(atoi(argv[3]));
	int rank_num;
	const int run_num(5);
	std::cout << "data_size:" << data_size << "  obj_num:" << obj_num /*<< "  rank_num:" << rank_num*/ << "  run_num:" << run_num << std::endl;
	//clock_t start(0), end(0), time_cost(0);
	std::chrono::time_point<std::chrono::system_clock> start_time;
	std::chrono::milliseconds time_cost;
	//NDS::new_uniform_population u1(obj_num, data_size, 0.5);
	//std::vector<std::vector<double>> data = u1.generate_output(rank_num);
	OFEC::uniform rand(0.5);
	NDS::circle_distribution u1(obj_num, data_size, rand);
	std::vector<std::vector<double>> data = u1.get_data();
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

	std::vector<int> YiyaSort_rank(data_size);
	std::vector<int> LinkSort_rank(data_size);
	std::vector<int> T_ENS_rank(data_size);
	std::vector<int> CornerSort_rank(data_size);
	std::vector<int> DeductiveSort_rank(data_size);
	std::vector<int> FastSort_rank(data_size);

	time_cost = time_cost.zero();
	for (int runID = 0; runID < run_num; ++runID) {
		start_time = std::chrono::system_clock::now();
		//start = clock();
		NDS::YiyaSort ys;
		ys.RankObjectiveLinkSort(data, YiyaSort_rank);
		//end = clock();
		//time_cost += end - start;
		time_cost += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);

	}
	std::cout << "\nYiyaSort:\t" << time_cost.count() / run_num << " milliseconds\n";

	time_cost = time_cost.zero();
	for (int runID = 0; runID < run_num; ++runID) {
		int ls_com = 0;
		start_time = std::chrono::system_clock::now();
		NDS::LinkSort(data, LinkSort_rank, ls_com);
		time_cost += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
	}
	std::cout << "\nLinkSort:\t" << time_cost.count() / run_num << " milliseconds\n";

	time_cost = time_cost.zero();
	for (int runID = 0; runID < run_num; ++runID) {
		int te_com = 0;
		start_time = std::chrono::system_clock::now();
		NDS::T_ENS(data, te_com, T_ENS_rank);
		time_cost += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
	}
	std::cout << "\nTree-ENS:\t" << time_cost.count() / run_num << " milliseconds\n";

	time_cost = time_cost.zero();
	for (int runID = 0; runID < run_num; ++runID) {
		double** POP = new double*[data_size];
		for (int i = 0; i < data_size; ++i)
			POP[i] = data[i].data();
		std::vector<int> cs_comp(data_size, 0);
		int num_comp;
		start_time = std::chrono::system_clock::now();
		NDS::cornerSort(POP, obj_num, data_size, CornerSort_rank.data(), cs_comp.data(), num_comp);
		time_cost += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
	}
	std::cout << "\nCornerSort:\t" << time_cost.count() / run_num << " milliseconds\n";

	time_cost = time_cost.zero();
	for (int runID = 0; runID < run_num; ++runID) {
		int ds_com = 0;
		start_time = std::chrono::system_clock::now();
		NDS::DeductiveSort(data, DeductiveSort_rank, ds_com);
		time_cost += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
	}
	std::cout << "\nDeductiveSort:\t" << time_cost.count() / run_num << " milliseconds\n";

	time_cost = time_cost.zero();
	for (int runID = 0; runID < run_num; ++runID) {
		start_time = std::chrono::system_clock::now();
		NDS::fast_sort fs(md);
		fs.sort();
		time_cost += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
		FastSort_rank = fs.rank_result();
		rank_num = fs.rank_num();
	}
	std::cout << "\nFastSort:\t" << time_cost.count() / run_num << " milliseconds\n" << std::endl;

	std::cout << "rank_num:" << rank_num << std::endl << std::endl;

	std::cout << "YiyaSort\t" << (FastSort_rank == YiyaSort_rank ? "succeed" : "failded") << std::endl;
	std::cout << "Tree_ENS\t" << (FastSort_rank == T_ENS_rank ? "succeed" : "failded") << std::endl;
	std::cout << "LinkSort\t" << (FastSort_rank == LinkSort_rank ? "succeed" : "failded") << std::endl;
	std::cout << "CornerSort\t" << (FastSort_rank == CornerSort_rank ? "succeed" : "failded") << std::endl;
	std::cout << "DeductiveSort\t" << (FastSort_rank == DeductiveSort_rank ? "succeed" : "failded") << std::endl;

	//system("pause");
	return 0;
}