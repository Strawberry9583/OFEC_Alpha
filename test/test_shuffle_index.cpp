#define BOOST_TEST_MODULE shuffle_index
#include <boost/test/unit_test.hpp>

#include<iostream>
#include "../utility/random/newran.h"

using namespace OFEC;

BOOST_AUTO_TEST_SUITE(shuffle_index_test)


BOOST_AUTO_TEST_CASE(test_case1) {
	uniform u(0.2);
	for (size_t i = 0; i < 10; i++) {
		std::cout << u.next() << "\t";
	}
	std::cout << std::endl;

	std::vector<double> v1(15);
	u.shuffle_index(v1, 12);
	for each (auto var in v1) {
		std::cout << var << "\t";
	}
	std::cout << std::endl;

	getchar();
}

BOOST_AUTO_TEST_SUITE_END()

