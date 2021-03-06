#pragma once
/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Li Zhou
* Email: 441837060@qq.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/

#ifndef N20_SR_RASTRIGIN_F9_H
#define N20_SR_RASTRIGIN_F9_H

#include "function_CEC2013.h"
namespace OFEC {
	namespace CEC2013 {
		class N20_SR_Rastrigin_F9 final:public function_CEC2013 {
		public:
			N20_SR_Rastrigin_F9(param_map &v);
			N20_SR_Rastrigin_F9(const std::string &name, size_t size_var, size_t size_obj);
			void evaluate__(real *x, std::vector<real>& obj);
			~N20_SR_Rastrigin_F9();
		protected:
			void initialize();
		};
	}
	using CEC2013_LSOP_F9 = CEC2013::N20_SR_Rastrigin_F9;
}
#endif


