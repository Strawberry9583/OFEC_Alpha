/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
TODO: compete the introduction of this class
*
*
*********************************************************************************/
#ifndef OFEC_REAL_DBG_H
#define OFEC_REAL_DBG_H

#include "../dynamic_continuous.h"
#include "../../../../../utility/matrix.h"

namespace OFEC {
	class real_DBG: public dynamic_continuous {
	protected:
		bool m_prediction;						// the next change of function can be predictable or not
		std::vector<matrix> mp_rotationMatrix;				// orthogonal rotation matrixes for each function
		std::vector<std::vector<std::vector<int>>> mppp_rotationPlanes;					// save the planes rotated during one periodicity

	public:
		real_DBG(const int size_var, const int number_peaks, const int size_object = 1);
		virtual ~real_DBG() = 0;

		///TODO: not 100% sure whether it should be derived by CompositionDBG and RotationDBG
		virtual bool set_period(const int period);

		real_DBG & operator=(const real_DBG & p);
		void reset();
		void reinitialize();
	protected:
		//TODO: delete this function?
		void correctSolution(double *x);
		void correctSolution(vector<double> x);
		void heightStandardChange();
		void positionStandardChange(double angle);

		virtual void randomChange() {};
		virtual void smallStepChange() {};
		virtual void largeStepChange() {};
		virtual void recurrentChange() {};
		virtual void chaoticChange() {};
		virtual void recurrentNoisyChange() {};


		void copy(problem * p);
		double  standardChange(const change_type T, const double min, const double max);
		virtual void allocateMemory(const int size_var, const int peaks);

		void initialize(const change_type t, const bool flag_var_change, const bool flag_number_peak_change,
			const int peak_number_change_mode, const bool flag_noise, const bool flag_time_linkage);
		void restoreInfor();
		virtual void changeDimension() {};
		virtual void changeNumPeaks() {};
	};
}

#endif // OFEC_REAL_DBG_H
