/*#define BOOST_TEST_MODULE CEC2013
#include <boost/test/unit_test.hpp>

#include "../instance/problem/continuous/large_scale/CEC2013/N7S1_ShiftedSchwefel_F7.h"
#include "../instance/problem/continuous/large_scale/CEC2013/ConflictingOS_ShiftedSchwefel_F14.h"
using namespace OFEC;

BOOST_AUTO_TEST_SUITE(CEC2013_test)


BOOST_AUTO_TEST_CASE(test_case1) {
	global::ms_global = std::unique_ptr<global>(new global(0.5, 0.5));
	std::vector<real> data(1000,0);
	std::vector<real> data2(1000,15);

	solution<variable<real>> sol(1, data);
	solution<variable<real>> sol2(1, data2);

	CEC2013::N7S1_ShiftedSchwefel_F7 a("N7S1_ShiftedSchwefel_F7", data.size(), 1);
	
	objective<real> temp_obj(1);
	solution<variable<real>, real> opt(a.get_optima().variable(0), std::move(temp_obj));
	a.evaluate(opt, caller::Problem);
	a.evaluate(sol, caller::Problem);
	a.evaluate(sol2, caller::Problem);
	for (auto &i : opt.get_objective())
		std::cout << i << std::endl;
	for (auto &i : sol.get_objective())
		std::cout << i << std::endl;
	for (auto &i : sol2.get_objective())
		std::cout << i << std::endl;

}
BOOST_AUTO_TEST_CASE(test_case2) {  // for F14
	global::ms_global = std::unique_ptr<global>(new global(0.5, 0.5));
	std::vector<real> data(1000, 0);
	std::vector<real> data2(1000, 15);

	solution<variable<real>> sol(1, data);
	solution<variable<real>> sol2(1, data2);

	CEC2013::ConflictingOS_ShiftedSchwefel_F14 a("ConflictingOS_ShiftedSchwefel_F14", data.size(), 1);

	objective<real> temp_obj(1);
	
	a.evaluate(sol, caller::Problem);
	a.evaluate(sol2, caller::Problem);
	
	std::cout << a.get_optima().single_objective(0) << std::endl;
	for (auto &i : sol.get_objective())
		std::cout << i << std::endl;
	for (auto &i : sol2.get_objective())
		std::cout << i << std::endl;

}



BOOST_AUTO_TEST_SUITE_END()*/