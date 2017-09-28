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
#ifndef OFEC_MOVINGPEAK_H
#define OFEC_MOVINGPEAK_H

#include "dynamic_continuous.h"

namespace OFEC {
	class moving_peak: public dynamic_continuous{
		int m_f;
		double m_vlength;					  /* distance by which the peaks are moved, severity */

		/* lambda determines whether there is a direction of the movement, or whether
		they are totally random. For lambda = 1.0 each move has the same direction,
		while for lambda = 0.0, each move has a random direction */
		double m_lambda;
		int m_use_basis_function;               /* if set to 1, a static landscape (basis_function) is included in the fitness evaluation */
		int m_calculate_right_peak;             /* saves computation time if not needed and set to 0 */
		double  m_standard_height;
		double  m_standard_width;              /* width chosen randomly when standardwidth = 0.0 */

		vector<double> m_shift;
		vector<int> m_covered_peaks;           /* which peaks are covered by the population ? */
		vector<vector<double>> m_prev_movement;/* to store every peak's previous movement */


	private:

		bool readData();

		// TODO: delete those functions with pointer paramater?
		/* the following basis functions are provided :*/
		double constant_basis_func(const double *gen);
		double constant_basis_func(const vector<double>& gen);
		double five_peak_basis_func(const double *gen);
		double five_peak_basis_func(const vector<double>& gen);
		/* the following peak functions are provided: */
		double peak_function1(const double *gen, int peak_number);
		double peak_function1(const std::vector<double>& gen, int peak_number);
		double peak_function_cone(const double *gen, const int &peak_number);
		double peak_function_cone(const std::vector<double>& gen, const int &peak_number);
		double peak_function_hilly(const double *gen, int peak_number);
		double peak_function_hilly(const std::vector<double>& gen, int peak_number);
		double peak_function_twin(const double  *gen, int peak_number);
		double peak_function_twin(const std::vector<double>&  gen, int peak_number);
		double function_selection(const double  *gen, const int &peak_number);
		double function_selection(const std::vector<double>& gen, const int & peak_number);
		double dummy_eval(const double *gen);
		void initialize();


	protected:
		void current_peak_calculate(const double *gen);
		virtual void allocate_memory(const int var_size, const int peaks);
		virtual void random_change();
		///TODO the flowing change types are not implemented
		virtual void small_step_change() { random_change(); };
		virtual void large_step_change() { random_change(); };
		virtual void recurrent_change() { random_change(); };
		virtual void chaotic_change() { random_change(); };
		virtual void recurrent_noisy_change() { random_change(); };
		virtual void change_dimension() { random_change(); };
		// TODO:kill pointer, and complete the pure virtual function
		virtual void change_num_peaks();
		//TODO: add "const" in parameter type.
		void copy(problem * problem);
		//16/05/2013
		virtual void update_peak_qaulity(); //TODO:???
		// TODO: should konw how to use the new evaluate()
		virtual void calculate_associate_radius();
		void set_severity();

	public:
		moving_peak(const int var_size, const int num_peaks, double const changing_ratio, const bool flag_var_change = false, \
			const bool flag_num_peak_change = false, const int peak_num_change_mode = 1, const bool flag_noise = false, const bool flag_time_linkage = false);
		moving_peak(param_map &v);
		//HACK: deleted it
		/*~moving_peak();*/

		/* assigns vlength a value from a normal distribution */
		void change_stepsize_random();
		/* sinusoidal change of the stepsize, */
		void change_stepsize_linear();
		/* returns 1 if current best individual is on highest m_peak, 0 otherwise */
		int get_right_peak();
		void set_v_length(const double s);
		void reset();
		moving_peak &operator=(moving_peak &other);
		//void reinitialize();
		double v_Length();
		//TODO: complete it
		// effective_fes: effective function evaluations constructed: is this problem object has been constructed
		evaluation_tag evaluate_(base &ss, caller call, bool effective_fes, bool constructed);
	};
}

#endif // !OFEC_MOVINGPEAK_H