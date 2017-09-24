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
* class dynamic_continuous defines dynamic continous optimization problems
*
*********************************************************************************/
#ifndef OFEC_DYNAMIC_CONTINUOUS_H
#define OFEC_DYNAMIC_CONTINUOUS_H

#include"dynamic.h"
#include"../continuous/continuous.h"

namespace OFEC {
	class dynamic_continuous : public dynamic, public continuous {
		
		// TODO: why don't store those mumbers of single peak to one class named peak? this class has so much data number now
		// TODO: What should I use, vector or unique_ptr?
		// int m_numPeaks;						        // number of peaks in Rotation_DBG , number of function in Composition_DBG
		vector<vector<double>> m_peak;          // positions of local or global optima(local optima in Rotation_DBG
		vector<vector<double>> m_pre_peak;		// global optima of basic function in Composition_DBG)
		vector<vector<double>> m_initial_peak;	// save the initial positions
		vector<double> m_height;					// peak height in Rotation_DBG, height of global optima in Composition_DBG
		vector<double> m_width;					    // weight value of each basic function in Composition_DBG,  peak width in Rotation_DBG

		///TODO preHeight and preWidth not considered in current version
		vector<double> m_pre_height;
		vector<double> m_pre_width;

		double m_min_height, m_max_height; // minimum\maximum height of all peaks(local optima) in Rotation_DBG(Composition_DBG)
		vector<double> m_height_severity;

		double m_min_width, m_max_width;
		vector<double> m_width_severity;

		vector<double> m_fit;						// objective value of each basic funciton in Composition_DBG, peak height in Rotation_DBG
		double m_global_optima;							// global optima value
		vector<bool> m_global_optima_idx;			// the index of the global optimal peak

		int m_current_peak;								// the peak where the best individual is located
		int m_max_peaks_number;							// the number of heigthest peaks
		double m_current_best;							// the objective value of current best individual

		vector<bool> m_whether_change;				// whether peaks change or not
		int m_num_change_peaks;							// the number of peaks that change
		double m_change_peak_ratio;						// the ratio of changing peaks

		int m_num_visable_peaks;                        // number of visable peaks, a peak is visable only if no peak is on top of it
		vector<int> m_is_tracked;					// accumulated number of peak[i] being tracked      
		vector<int> m_height_order;
		int m_peaks_found;
		vector<bool> m_found;
		vector<double> m_time_linkage;				// a random vector to change peaks position which are being tracked
		vector<int> m_amended_height_order;
		//added 16/05/2013
		vector<double> m_associate_radius; /*actual radius in the fitness landscape*/
		double m_peak_qaulity;	/*to evaluate qaulity of peaks trcked  in terms of peaks heights*/
		//added 04/07/2014
		bool is_global_opt(int idx); // TODO: is global optima exists in dynamic problem?
	public:
		dynamic_continuous(const int size_var, const int num_peaks, const unsigned size_obj = 1);
		virtual ~dynamic_continuous() = 0;

		double get_global_max()const;
		void print_fun(std::ofstream & out);
		const vector<double>& get_peak(const int p) const;
		const vector<vector<double>>& get_all_peaks()const; // TDOD:???
		double get_peak_height(const int p)const;
		double get_pre_peak_height(const int p)const;
		double get_pre_peak_width(const int p)const;
		const vector<double>& get_height() const;
		const vector<double>& get_pre_peak(const int p)const;
		const vector<bool>& get_global_optima_idx()const;
		int get_number_of_global_opt_peak()const;

		void set_number_of_changes(const int n);
		void set_number_of_changes(const double ratio);

		void set_height_severity(const double s);
		void set_width_severity(const double s);
		dynamic_continuous &operator=(const dynamic_continuous &dynamic_continuous_pro);

		void set_height(const vector<double>& h);
		void set_position(const vector<vector<double>>& p);//const double **p 
		virtual void set_width(const double w);

		/// for debug mode
		void print_peak(const int idx);
		void print_peaks(ofstream & out);
		int get_num_of_visable_peaks();
		bool is_visable(const int idx); // TODO: Should I take out the "is_" in the function name?
		int get_track_number(int idex);
		bool is_tracked(vector<double> &gen, vector<double> &obj); // is any peak tracked for the first time
		bool is_tracked(double *gen, vector<double> &obj); //TODO:Should I use template to recompose the two functions
		int get_peaks_found();
		double get_associate_radius(int idx);
		double get_peaks_traced_qaulity();
		//15-07-2013
		bool is_global_optima_tracked();
		// TODO: Need a class my_vector here
		const double * get_nearest_peak(const vector<double>& p);
	protected:
		virtual void random_change() {};
		virtual void small_step_change() {};
		virtual void large_step_change() {};
		virtual void recurrent_change() {};
		virtual void chaotic_change() {};
		virtual void recurrent_noisy_change() {};

		void copy(problem * p); // copy parameter values of a problem when it changes
		//HACK: deleted it
		//virtual void free_memory();
		virtual void allocate_memory(const int size_var, const int num_peaks);

		virtual void change_variable() {};
		virtual void change_num_peaks() {};

		void calculate_global_optima();
		void update_number_of_changes();
		void compute_num_visable_peaks();
		void add_noise(double *x_);
		void update_time_linkage();
		void move_peak(const int idx);
		virtual void calculate_associate_radius() {};
		virtual void update_peak_qaulity() {};
	};
}

#endif // !OFEC_DYNAMIC_CONTINUOUS_H