/*#define BOOST_TEST_MODULE CEC2015_expensive
#include <boost/test/unit_test.hpp>

#include "../instance/problem/continuous/expensive/CEC2015/F6_shifted_rotated_happy_cat.h"

#include "../instance/problem/continuous/expensive/CEC2015/F10_hybrid1.h"
#include "../instance/problem/continuous/expensive/CEC2015/F13_expensive_composition1.h"
using namespace OFEC;

BOOST_AUTO_TEST_SUITE(CEC2015_expensive_test)


BOOST_AUTO_TEST_CASE(test_case1) {
	global::ms_global = std::unique_ptr<global>(new global(0.5, 0.5));

	CEC2015::F6_shifted_rotated_happy_cat a("F6_shifted_rotated_happy_cat", 10, 1);

	//vector<real> data = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	std::vector<real> data(10, 0);
	if(a.translation().size()>0)
		for (size_t i = 0; i < data.size(); ++i) {
			data[i] += a.translation()[i];
		}
	
	solution<variable<real>> sol(1, data);

	a.evaluate(sol, caller::Problem);
	BOOST_CHECK_EQUAL(a.get_optima().single_objective(0), sol.get_objective()[0]);

	std::cout << a.get_optima().single_objective(0) << std::endl;

}

BOOST_AUTO_TEST_CASE(test_case2) {
	global::ms_global = std::unique_ptr<global>(new global(0.5, 0.5));

	CEC2015::F10_hybrid1 a("F10_hybrid1", 10, 1);

	//vector<real> data = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	std::vector<real> data(10, 0);

	
	//for (size_t i = 0; i < data.size(); ++i) {
	//	data[i] = a.get_function(0)->get_optima().variable(0)[i];
	//}
	solution<variable<real>> sol(1, data);

	a.evaluate(sol, caller::Problem);
	BOOST_CHECK_EQUAL(a.get_optima().single_objective(0), sol.get_objective()[0]);

	std::cout << a.get_optima().single_objective(0) << std::endl;
	
}
BOOST_AUTO_TEST_CASE(test_case3) {
	global::ms_global = std::unique_ptr<global>(new global(0.5, 0.5));

	CEC2015::F13_expensive_composition1 a("F13_expensive_composition1", 10, 1);

	//vector<real> data = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	std::vector<real> data(10, 1);

	
	for (size_t i = 0; i < data.size(); ++i) {
		data[i] += a.get_function(0)->translation()[i];
	}
	solution<variable<real>> sol(1, data);

	a.evaluate(sol, caller::Problem);
	BOOST_CHECK_EQUAL(a.get_optima().single_objective(0), sol.get_objective()[0]);

	std::cout << a.get_optima().single_objective(0) << std::endl;

}

BOOST_AUTO_TEST_SUITE_END()*/