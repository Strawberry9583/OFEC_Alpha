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

#include "scaffer_F6.h"
namespace OFEC {
	
	scaffer_F6::scaffer_F6(param_map &v) :problem((v[param_proName]), (v[param_numDim]), 1), \
		function((v[param_proName]), (v[param_numDim]), 1) {

		set_range(-100, 100);
		set_init_range(-100, 100);
		initialize();
	}
	scaffer_F6::scaffer_F6(const std::string &name, size_t size_var, size_t size_obj) :problem(name, size_var, size_obj), \
		function(name, size_var, size_obj) {

		set_range(-100, 100);
		set_init_range(-100, 100);
		initialize();
	}

	void scaffer_F6::initialize() {

		set_original_global_opt();
		set_bias(0);
		set_global_opt();

		m_variable_accuracy = 1.0e-2;
	}

	void scaffer_F6::evaluate__(real *x, std::vector<real>& obj) {
		if (m_translation_flag)
			translate(x);
		if (m_scale_flag)
			scale(x);
		if (m_rotation_flag)
			rotate(x);
		if (m_translation_flag)
			translate_origin(x);
		double fitness = 0;
		for (size_t i = 0; i < m_variable_size-1; ++i) {
			fitness += 0.5 + (sin(sqrt(x[i] * x[i] + x[i+1] * x[i+1]))*sin(sqrt(x[i] * x[i] + x[i+1] * x[i+1])) - 0.5) / ((1 + 0.001*(x[i] * x[i] + x[i+1] * x[i+1]))*(1 + 0.001*(x[i] * x[i] + x[i+1] * x[i+1])));
		}
		obj[0] = fitness + m_bias;

	}
	
}