/*boost unit test for multi-objective non-dominated sorting*/

#define BOOST_TEST_MODULE non_dominated_sorting
#include <boost/test/unit_test.hpp>
#include "../utility/nondominated_sorting/DominanceGraphSort.h"
#include "../utility/nondominated_sorting/FastNondominatedSort.h"
#include "../utility/nondominated_sorting/CornerSort.h"
#include "../utility/nondominated_sorting/preSorting.h"
#include "../utility/nondominated_sorting/nodes_initialize.h"
#include "../utility/nondominated_sorting/T_ENS.h"
#include "../utility/nondominated_sorting/LinkSort.h"
#include <time.h>
#include <fstream>

BOOST_AUTO_TEST_SUITE(nondominated_sorting_test)

/*
BOOST_AUTO_TEST_CASE(test_case4) {
	const int data_size(1000);
	const int obj_num(5);
	const int test_number(50);
	const double change_ratio(0.5);
	random rand(0.5);
	ofstream file("e:\\50percent.dat");

	for (int ith = 0; ith < test_number; ++ith) {
		cout << ith + 1 << endl;
		NDS::circle_distribution distribution(obj_num, data_size, rand);
		vector<vector<double>> data(distribution.get_data());

		//FNS
		vector<int> no;
		vector<pair<int, vector<double>>> md;
		for (int i = 0; i < data_size; ++i) {
			no.push_back(i);
			md.emplace_back(no[i], data[i]);
		}

		NDS::fast_sort fs(md);
		fs.sort();
		cout << fs.number() << "\t";
		file << fs.number() << "\t";

		//CS
		double** POP = new double*[data_size];
		for (int i = 0; i < data_size; ++i) {
			POP[i] = new double[obj_num];
			for (int j = 0; j < obj_num; ++j)
				POP[i][j] = data[i][j];
		}
		unsigned int cs_rank[data_size];
		int cs_num(0);
		int cs_com[obj_num] = { 0 };
		NDS::cornerSort(POP, obj_num, data_size, cs_rank, cs_com, cs_num);
		cout << cs_num << "\t";
		file << cs_num << "\t";

		//DG(Rebuild)
		NDS::dg_graph<double> dg(data);
		cout << dg.number() << "\t";
		file << dg.number() << "\t";

		//DG(Replace)
		cout << endl;
		file << endl;

		//Change
		int temp_DG_Replace_number = dg.number();

		for (int j = data_size - 1; j > int(data_size * (1 - change_ratio)) - 1; --j) {
			dg._delete(j);
			data.pop_back();
		}
		for (int j = data_size - 1; j > int(data_size * (1 - change_ratio)) - 1; --j) {
			auto temp_new_one = distribution.new_one(rand);
			data.push_back(temp_new_one);
			dg.insert(&data.back());
		}

		//FNS
		no.clear();
		md.clear();
		for (int i = 0; i < data_size; ++i) {
			no.push_back(i);
			md.emplace_back(no[i], data[i]);
		}
		NDS::fast_sort fs_c(md);
		fs_c.sort();
		cout << fs_c.number() << "\t";
		file << fs_c.number() << "\t";

		//CS
		double** POP_c = new double*[data_size];
		for (int i = 0; i < data_size; ++i) {
			POP_c[i] = new double[obj_num];
			for (int j = 0; j < obj_num; ++j)
				POP_c[i][j] = data[i][j];
		}
		unsigned int cs_rank_c[data_size];
		int cs_num_c(0);
		NDS::cornerSort(POP_c, obj_num, data_size, cs_rank_c, cs_com, cs_num_c);
		cout << cs_num_c << "\t";
		file << cs_num_c << "\t";

		//DG(Rebuild)
		NDS::dg_graph<double> dg_c(data);
		cout << dg_c.number() << "\t";
		file << dg_c.number() << "\t";

		//DG(Replace)
		cout << dg.number() - temp_DG_Replace_number << endl;
		file << dg.number() - temp_DG_Replace_number << endl;
	}
	file.close();
	system("pause");

}
*/


/*
BOOST_AUTO_TEST_CASE(test_case3) {
	const int data_size(1000);
	const int obj_num(3);

	NDS::circle_distribution d1(obj_num, data_size);
	vector<vector<double>> data = d1.get_data();
	vector<int> no;
	vector<pair<int, vector<double>>> md;
	for (int i = 0; i < data_size; ++i) {
		no.push_back(i);
		md.emplace_back(no[i], data[i]);
	}
	NDS::d_quicktree dqt(obj_num);
	dqt.build(md);
	dqt._delete(149);
	dqt.add(vector<double>({ 0.1,0.1,0.1 }));
	auto dqt_result = dqt.ranking();

	data.erase(data.begin() + 149);
	data.push_back(vector<double>({ 0.1,0.1,0.1 }));
	no.clear();
	md.clear();
	for (int i = 0; i < data_size; ++i) {
		no.push_back(i);
		md.emplace_back(no[i], data[i]);
	}
	NDS::fast_sort fs(md);
	fs.sort();
	auto & fs_result = fs.ranking();

	for (auto x : fs_result)
		if (x.second == 0) cout << x.first << " ";	cout << endl;

	for (auto x : dqt_result)
		if (x.second == 0) cout << x.first << " ";	cout << endl;

	system("pause");
}
*/
BOOST_AUTO_TEST_CASE(test_case2) {
	const int data_size(1000);
	const int obj_num(3);
	OFEC::random rand1(0.5);

	NDS::circle_distribution d1(obj_num, data_size, rand1);
	std::vector<std::vector<double>> data = d1.get_data();
	std::vector<int> no;
	std::vector<std::pair<int, std::vector<double>>> md;
	for (int i = 0; i < data_size; ++i) {
		no.push_back(i);
		md.emplace_back(no[i], data[i]);
	}

	NDS::fast_sort fs(md);
	fs.sort();
	auto fs_result = fs.ranking();
	int fs_numcom = fs.number();

	std::vector<std::vector<double>> data2(obj_num, std::vector<double>(data_size));
	for (int i = 0; i < data_size; ++i)
		for (int j = 0; j < obj_num; ++j)
			data2[j][i] = data[i][j];

	int numcom = 0;
	std::map<int, int> ls_result;
	NDS::LinkSort(data2, ls_result, numcom);

	std::cout << (ls_result == fs_result ? "match" : "not match") << std::endl;
	system("pause");
	
}
/*
BOOST_AUTO_TEST_CASE(test_case1) {

	const int data_size(1000);
	const int obj_num(5);

	time_t start(0), end(0);
	NDS::line_distribution d1(obj_num, data_size);
	vector<vector<double>> data = d1.get_data();
	vector<int> no;
	vector<pair<int, vector<double>>> md;
	double** POP = new double*[data_size];
	for (int i = 0; i < data_size; ++i) {
		no.push_back(i);
		md.emplace_back(no[i], data[i]);
		POP[i] = new double[obj_num];
		for (int j = 0; j < obj_num; ++j)
			POP[i][j] = data[i][j];
	}
	unsigned int cs_rank[data_size]; 
	int cs_num(0);
	int cs_com[obj_num] = {0};
	map<int, int> cs_result;


	cout << "NSGA¢ò Non-dominated sorting:" << endl;
	time(&start);
	NDS::fast_sort fs(md);
	fs.sort();
	auto & fs_result = fs.ranking();
	time(&end);
	cout << "Compares " << fs.number() << " times" << endl;
	cout << "cost " << (end - start) << " seconds" << endl << endl;

	cout << "Corner Sort Non-dominated sorting:" << endl;
	time(&start);
	NDS::cornerSort(POP, obj_num, data_size, cs_rank, cs_com, cs_num);
	for (int i = 0; i < data_size; ++i)
		cs_result.emplace(i, cs_rank[i]-1);
	time(&end);
	cout << "Compares " << cs_num << " times" << endl;
	cout << "cost " << (end - start) << " seconds" << endl << endl;
	
	cout << "Quick Dominance Tree Non-dominated sorting:" << endl;
	time(&start);
	NDS::d_quicktree dqt(obj_num);
	dqt.build(md);
	auto dqt_result = dqt.ranking();
	time(&end);
	cout << "Compares " << dqt.number() << " times" << endl;
	cout << "cost " << (end - start) << " seconds" << endl << endl;

	cout << "Dominance graph Non-dominated sorting:" << endl;
	time(&start);
	NDS::dg_graph<double> dg(data);
	auto dg_result = dg.ranking();
	time(&end);
	cout << "Compares " << dg.number() << " times" << endl;
	cout << "Cost " << (end - start) << " seconds" << endl << endl;
	
	cout << ((fs_result == cs_result && cs_result ==  dg_result) ? "match" : "not match") << endl;
	system("pause");
}
*/

BOOST_AUTO_TEST_SUITE_END()