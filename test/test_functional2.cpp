#define BOOST_TEST_MODULE functional
#include <boost/test/unit_test.hpp>

#include<iostream>
#include "../utility/random/newran.h"
#include"../utility/functional.h"

using namespace OFEC;

// test the shuffle_index() and get_random() in functional.h
BOOST_AUTO_TEST_SUITE(functional_test)


BOOST_AUTO_TEST_CASE(test_get_random) {

	uniform u(0.2);

	std::cout << "mean: " << u.mean() << std::endl;
	for (size_t i = 0; i < 10; i++) {
		std::cout << u.next() << "\t";
	}
	std::cout << std::endl;

	std::cout << "get float random: " << std::endl;
	for (size_t i = 0; i < 10; i++) {
		std::cout<<get_random(0.0, 10.0, &u)<<"\t";
	}
	std::cout << std::endl;

	std::cout << "get int random" << std::endl;
	for (size_t i = 0; i < 10; i++) {
		std::cout << get_random(0, 10, &u) << "\t";
	}
	std::cout << std::endl << std::endl;
}

BOOST_AUTO_TEST_CASE(test_shuffle_index) {

	std::vector<int> v1(15);
	uniform u(0.2);
	unique_ptr<uniform> u1 = make_unique<uniform>(u);
	shuffle_index(v1, 10, u1.get());
	std::cout << "the shuffled index: " << std::endl;
	for each (int var in v1) {
		std::cout << var << "\t";
	}
	std::cout << std::endl;

	getchar();
}

BOOST_AUTO_TEST_SUITE_END()


