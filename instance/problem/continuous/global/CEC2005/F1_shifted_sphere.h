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


#ifndef OFEC_F1_SHIFTED_SPHERE_H
#define OFEC_F1_SHIFTED_SPHERE_H

#include "../classical/sphere.h"

namespace OFEC {
	namespace CEC2005 {
		class F1_shifted_sphere final : public sphere
		{
		public:
			F1_shifted_sphere(param_map &v);
			F1_shifted_sphere(const std::string &name, size_t size_var, size_t size_obj);
		protected:
			void initialize();
			void evaluate__(real *x, std::vector<real>& obj);
		private:
		};
	}
	using CEC2005_GOP_F1 = CEC2005::F1_shifted_sphere;
}
#endif // ! OFEC_F1_SHIFTED_SPHERE_H
