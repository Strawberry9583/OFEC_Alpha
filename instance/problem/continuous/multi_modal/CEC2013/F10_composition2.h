/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li and Li Zhou
* Email: changhe.lw@gmail.com, 441837060@qq.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************
*  Paper: Benchmark Functions for CEC��2013 Special Session and Competition on Niching
*  Methods for Multimodal Function Optimization.
*******************************************************************************************/

#ifndef OFEC_F10_COMPOSITION2_H
#define OFEC_F10_COMPOSITION2_H

#include "../../global/CEC2005/composition.h"

namespace OFEC {
	namespace CEC2013 {
		class F10_composition2 final : public CEC2005::composition
		{
		public:
			F10_composition2(param_map &v);
			F10_composition2(const std::string &name, size_t size_var, size_t size_obj);
		protected:
			void initialize();
			void evaluate__(real *x, std::vector<real>& obj);
			void set_function();

			void set_rotation();
		private:

		};
	}
	using CEC2013_MMP_F10 = CEC2013::F10_composition2;
}
#endif // !OFEC_F10_COMPOSITION2_H


