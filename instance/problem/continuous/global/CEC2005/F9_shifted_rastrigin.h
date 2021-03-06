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
*************************************************************************/


#ifndef OFEC_F9_SHIFTED_RASTRIGIN_H
#define OFEC_F9_SHIFTED_RASTRIGIN_H

#include "../classical/rastrigin.h"

namespace OFEC {
	namespace CEC2005 {
		class F9_shifted_rastrigin final: public rastrigin
		{
		public:
			F9_shifted_rastrigin(param_map &v);
			F9_shifted_rastrigin(const std::string &name, size_t size_var, size_t size_obj);
		protected:
			void initialize();
			void evaluate__(real *x, std::vector<real>& obj);
		private:
		};
	}
	using CEC2005_GOP_F9 = CEC2005::F9_shifted_rastrigin;
}
#endif // ! OFEC_F9_SHIFTED_RASTRIGIN_H


