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
#ifndef OFEC_COMPOSITION_DBG
#define OFEC_COMPOSITION_DBG

#include"../real_DBG.h"

namespace OFEC {
	enum class function_tag { Sphere = 0, Rastrigin, Griewank, Ackley, Weierstrass };
	enum class composition_DBG_function_id { COMDBG_SPHERE = 0, COMDBG_RASTRIGIN, COMDBG_GRIEWANK, COMDBG_ACKLEY, COMDBG_HYBRID };
	class composition_DBG :public real_DBG {
	private:
		//TODO: the original Boundary looks like domain<real>
		vector<domain<real>> m_com_boundary;				    // boundary of component functions
		std::vector<double> m_converge_severity;				// severity of converge range for each function
		std::vector<double> m_stretch_severity;					// severity of stretching original function, greater than 1 for stretch
																// less than 1 for compress the original function
		double m_height_normalize_severity;					    // constant number for noralizing all basic function with similar height
		std::vector<function_tag> m_com_function;						        // which basic function used to compose the composition function

		std::vector<double> m_fmax;
		static const int m_num_com_funs = 5;					// number of basic functions
		static thread_local unique_ptr<composition_DBG_function_id> ms_fun_id;
	private:

		void set_com_boundary();
		void get_com_boundary(const function_tag &f, double &l, double &u, const int rDimIdx = 0)const;
		void initialize(const change_type T, const composition_DBG_function_id rF, const double rChangingRatio, const bool rFlagDimChange = false, \
			const bool rFlagNumPeakChange = false, const int peakNumChangeMode = 1, const bool flagNoise = false, const bool flagTimelinkage = false);

	public:
		static const int msc_num_funs = 5;                     // number of functions in comDBG system
		composition_DBG(const int rDimNumber, const int rNumPeaks, const change_type rT, const composition_DBG_function_id rF, \
			const double rChangingRatio, const bool rFlagDimChange, const bool rFlagNumPeakChange, \
			const int peakNumChangeMode, const bool flagNoise, const bool flagTimelinkage);
		composition_DBG(param_map &v);

		virtual ~composition_DBG();

		composition_DBG &operator =(const composition_DBG &);
		void set_rotation_matrix();					        //randomly generate rotation matrix for each basic function
		void set_coverge_sevrity(const double* cs);
		void set_stretch_severity();
		void set_basic_function(const function_tag *bf);		//component functions to compose the search space	
		evaluation_tag evaluate_(base &ss, caller call, bool effective_eval, bool constructed);
	protected:
		void copy(const problem * rP);
		virtual void  free_memory();

		virtual void random_change();
		virtual void small_step_change();
		virtual void large_step_change();
		virtual void recurrent_change();
		virtual void chaotic_change();
		virtual void recurrent_noisy_change();
		virtual void allocate_memory(const int rDimNum, const int rPeaks);
		virtual void change_variable();
		virtual void change_num_peaks();

	private:
		void correct_solution(const function_tag &f, double *x);	// make x within search range after rotation											// basic five functions
		double f_sphere(double *x);
		double f_rastrigin(double *x);
		double f_weierstrass(double *x);
		double f_griewank(double *x);
		double f_ackley(double *x);
		double select_fun(const function_tag &f, double *x);
	};

}

#endif // !OFEC_COMPOSITION_DBG

