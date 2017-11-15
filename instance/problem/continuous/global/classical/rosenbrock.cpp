/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com 
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/

#include "rosenbrock.h"

namespace OFEC {
	
	rosenbrock::rosenbrock(param_map &v) :problem((v[param_proName]), (v[param_numDim]), 1), \
		function((v[param_proName]), (v[param_numDim]), 1) {

		set_range(-2.048, 2.048);
		set_init_range(-2.048, 2.048);
		initialize();
	}
	rosenbrock::rosenbrock(const std::string &name, size_t size_var, size_t size_obj) :problem(name, size_var, size_obj), \
		function(name, size_var, size_obj) {

		set_range(-2.048, 2.048);
		set_init_range(-2.048, 2.048);
		initialize();
	}

	void rosenbrock::initialize() {
		std::vector<double> v(m_variable_size, 1);
		set_original_global_opt(v.data());

		set_global_opt();
		//m_objective_accuracy = 1.0e-5;
	}

	void rosenbrock::evaluate__(real *x, std::vector<real>& obj) {
		if (m_translation_flag)
			translate(x);
		if (m_scale_flag)
			scale(x);
		if (m_rotation_flag)
			rotate(x);
		if (m_translation_flag)
			translate_origin(x);
		double fitness = 0;

		for (int i = 0; i < m_variable_size - 1; i++) {
			fitness += 100 * (x[i + 1] - x[i] * x[i])*(x[i + 1] - x[i] * x[i]) + (x[i] - 1)*(x[i] - 1);
		}

		obj[0] = fitness + m_bias;

	}
	
}