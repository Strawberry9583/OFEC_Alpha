#include "FastSort.h"
#include "CornerSort.h"
#include "T_ENS.h"
#include "LinkSort.h"
#include "YiyaSort.h"
#include "DeductiveSort.h"
#include "new_static_population.h"
#include "nodes_initialize.h"

#include <chrono>
#include <fstream>

int main(int argc, char* argv[]) {
	int data_size(atoi(argv[1]));
	int obj_num(atoi(argv[2]));
	//int rank_num(atoi(argv[3]));
	int rank_num;
	const int run_num(5);
	std::cout << "data_size:" << data_size << "  obj_num:" << obj_num /*<< "  rank_num:" << rank_num*/ << "  run_num:" << run_num << std::endl;
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

	std::vector<int> YiyaSort_rank(data_size);
	std::vector<int> LinkSort_rank(data_size);
	std::vector<int> T_ENS_rank(data_size);
	std::vector<int> CornerSort_rank(data_size);
	std::vector<int> DeductiveSort_rank(data_size);
	std::vector<int> FastSort_rank(data_size);

	time_cost = time_cost.zero();
	for (int runID = 0; runID < run_num; ++runID) {
		start_time = std::chrono::system_clock::now();
		NDS::YiyaSort ys;
		ys.RankObjectiveLinkSort(data, YiyaSort_rank);
		time_cost += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
	}
	std::cout << "YiyaSort:\t" << time_cost.count() / run_num << " ms" << std::endl;

	time_cost = time_cost.zero();
	for (int runID = 0; runID < run_num; ++runID) {
		int ls_com = 0;
		start_time = std::chrono::system_clock::now();
		NDS::LinkSort(data, LinkSort_rank, ls_com);
		time_cost += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
	}
    #ifdef USING_CONCURRENT
	std::cout << "LinkSort(CON):\t" << time_cost.count() / run_num << " ms" << std::endl;
    #else
	std::cout << "LinkSort:\t" << time_cost.count() / run_num << " ms" << std::endl;
    #endif // USING_CONCURRENT

	time_cost = time_cost.zero();
	for (int runID = 0; runID < run_num; ++runID) {
		int te_com = 0;
		start_time = std::chrono::system_clock::now();
		NDS::T_ENS(data, te_com, T_ENS_rank);
		time_cost += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
	}
	std::cout << "Tree-ENS:\t" << time_cost.count() / run_num << " ms" << std::endl;

	time_cost = time_cost.zero();
	for (int runID = 0; runID < run_num; ++runID) {
		start_time = std::chrono::system_clock::now();
		int num_comp(0);
		NDS::CornerSort(data, CornerSort_rank, num_comp);
		time_cost += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
	}
	std::cout << "CornerSort:\t" << time_cost.count() / run_num << " ms" << std::endl;

	time_cost = time_cost.zero();
	for (int runID = 0; runID < run_num; ++runID) {
		int ds_com = 0;
		start_time = std::chrono::system_clock::now();
		NDS::DeductiveSort(data, DeductiveSort_rank, ds_com);
		time_cost += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
	}
	std::cout << "DeductiveSort:\t" << time_cost.count() / run_num << " ms" << std::endl;

	time_cost = time_cost.zero();
	for (int runID = 0; runID < run_num; ++runID) {
		int num_comp(0);
		start_time = std::chrono::system_clock::now();
		rank_num = NDS::FastSort(data, FastSort_rank, num_comp);
		time_cost += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
	}
    #ifdef USING_CONCURRENT
	std::cout << "FastSort(CON):\t" << time_cost.count() / run_num << " ms" << std::endl;
    #else
	std::cout << "FastSort:\t" << time_cost.count() / run_num << " ms" << std::endl;
    #endif // USING_CONCURRENT

	std::cout << "rank_num:" << rank_num << std::endl;

	std::cout << "YiyaSort\t" << (FastSort_rank == YiyaSort_rank ? "correct" : "incorrect") << std::endl;
	std::cout << "Tree_ENS\t" << (FastSort_rank == T_ENS_rank ? "correct" : "incorrect") << std::endl;
	std::cout << "LinkSort\t" << (FastSort_rank == LinkSort_rank ? "correct" : "incorrect") << std::endl;
	std::cout << "CornerSort\t" << (FastSort_rank == CornerSort_rank ? "correct" : "incorrect") << std::endl;
	std::cout << "DeductiveSort\t" << (FastSort_rank == DeductiveSort_rank ? "correct" : "incorrect") << std::endl;

	return 0;
}
