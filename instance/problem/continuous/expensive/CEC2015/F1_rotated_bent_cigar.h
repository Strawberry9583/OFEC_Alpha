
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


#ifndef OFEC_F1_ROTATED_BENT_CIGAR_H
#define OFEC_F1_ROTATED_BENT_CIGAR_H

#include "../../global/classical/bent_cigar.h"

namespace OFEC {
	namespace CEC2015 {
		class F1_rotated_bent_cigar final : public bent_cigar
		{
		public:
			F1_rotated_bent_cigar(param_map &v);
			F1_rotated_bent_cigar(const std::string &name, size_t size_var, size_t size_obj);
		protected:
			void initialize();
			void evaluate__(real *x, std::vector<real>& obj);
		private:
		};
	}
	using CEC2015_EOP_F1 = CEC2015::F1_rotated_bent_cigar;
}
#endif // ! OFEC_F1_ROTATED_BENT_CIGAR_H
