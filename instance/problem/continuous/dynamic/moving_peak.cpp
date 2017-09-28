#include "moving_peak.h"

namespace OFEC {

	static double basis_peak[5][7] = {
		{ 8.0,  64.0,  67.0,  55.0,   4.0, 0.1, 50.0 },
		{ 50.0,  13.0,  76.0,  15.0,   7.0, 0.1, 50.0 },
		{ 9.0,  19.0,  27.0,  67.0, 24.0, 0.1, 50.0 },
		{ 66.0,  87.0,  65.0,  19.0,  43.0, 0.1, 50.0 },
		{ 76.0,  32.0,  43.0,  54.0,  65.0, 0.1, 50.0 }
	};

	static double twin_peak[7] = /* difference to first peak */
	{
		1.0,  1.0,  1.0,  1.0,   1.0, 0.0, 0.0,
	};


	bool moving_peak::readData() {
		/***************************************
		//		m_f		Evaluation Function
		//		1		constant_basis_func()
		//		2		five_peak_basis_func()
		//		3		peak_function1()
		//		4		peak_function_cone()
		//		5		peak_function_hilly()
		//		6		peak_function_twin()
		**************************************
		in>>temp>>m_vlength; // distance by which the peaks are moved, severity
		lambda determines whether there is a direction of the movement, or whether
		they are totally random. For lambda = 1.0 each move has the same direction,
		while for lambda = 0.0, each move has a random direction
		//in>>temp>>lambda;
		//in>>temp>>m_use_basis_function;  if set to 1, a static landscape (basis_function) is included in the fitness evaluation
		}*/

		m_f = 4;
		m_lambda = 0;
		m_use_basis_function = 0;
		m_vlength = 1.0;
		m_calculate_right_peak = 1; /* saves computation time if not needed and set to 0 */
		m_min_height = 30.0;
		m_max_height = 70.0;
		m_standard_height = 50.0;
		m_min_width = 1.0;
		m_max_width = 12.0;
		m_standard_width = 0.0;

		set_severity();
		return true;
	}

	double moving_peak::constant_basis_func(const double * gen) {
		return 0.0;
	}
	double moving_peak::constant_basis_func(const vector<double>& gen) {
		return 0.0;
	}

	//TODO: delete it?
	double moving_peak::five_peak_basis_func(const double * gen) {
		double maximum = -100000.0, dummy = 0;
		for (int i = 0; i < 5; i++) {
			dummy = (gen[0] - basis_peak[i][0])*(gen[0] - basis_peak[i][0]);
			for (int j = 1; j < m_variable_size; j++)  dummy += (gen[j] - basis_peak[i][j])*(gen[j] - basis_peak[i][j]);
			dummy = basis_peak[i][m_variable_size + 1] - (basis_peak[i][m_variable_size] * dummy);
			if (dummy > maximum)       maximum = dummy;
		}
		return maximum;
	}
	double moving_peak::five_peak_basis_func(const vector<double>& gen) {
		double maximum = -100000.0, dummy = 0;
		for (int i = 0; i < 5; i++) {
			dummy = (gen[0] - basis_peak[i][0])*(gen[0] - basis_peak[i][0]);
			for (int j = 1; j < m_variable_size; j++)  dummy += (gen[j] - basis_peak[i][j])*(gen[j] - basis_peak[i][j]);
			dummy = basis_peak[i][m_variable_size + 1] - (basis_peak[i][m_variable_size] * dummy);
			if (dummy > maximum)       maximum = dummy;
		}
		return maximum;
	}

	//TODO: can I delete this function?
	double moving_peak::peak_function1(const double * gen, int peak_number) {
		double dummy = (gen[0] - m_peak[peak_number][0])*(gen[0] - m_peak[peak_number][0]);
		for (int j = 1; j < m_variable_size; j++)
			dummy += (gen[j] - m_peak[peak_number][j])*(gen[j] - m_peak[peak_number][j]);
		return m_height[peak_number] / (1 + m_width[peak_number] * dummy);
	}
	double moving_peak::peak_function1(const std::vector<double>& gen, int peak_number) {
		double dummy = (gen[0] - m_peak[peak_number][0])*(gen[0] - m_peak[peak_number][0]);
		for (int j = 1; j < m_variable_size; j++)
			dummy += (gen[j] - m_peak[peak_number][j])*(gen[j] - m_peak[peak_number][j]);
		return m_height[peak_number] / (1 + m_width[peak_number] * dummy);
	}

	//TODO: delete it?
	double moving_peak::peak_function_cone(const double * gen, const int & peak_number) {
		double val, dummy = 0;
		for (int j = 0; j < m_variable_size; j++) {
			val = gen[j] - m_peak[peak_number][j];
			dummy += val*val;
		}
		if (dummy != 0)  dummy = m_height[peak_number] - m_width[peak_number] * sqrt(dummy);
		return dummy;
	}
	double moving_peak::peak_function_cone(const std::vector<double>& gen, const int & peak_number) {
		double val, dummy = 0;
		for (int j = 0; j < m_variable_size; j++) {
			val = gen[j] - m_peak[peak_number][j];
			dummy += val*val;
		}
		if (dummy != 0)  dummy = m_height[peak_number] - m_width[peak_number] * sqrt(dummy);
		return dummy;
	}

	double moving_peak::peak_function_hilly(const double * gen, int peak_number) {
		int j = 0;
		double dummy = (gen[0] - m_peak[peak_number][0])*(gen[0] - m_peak[peak_number][0]);
		for (j = 1; j < m_variable_size; j++)
			dummy += (gen[j] - m_peak[peak_number][j])*(gen[j] - m_peak[peak_number][j]);

		return m_height[peak_number] - m_width[peak_number] * dummy - 0.01*sin(20.0*dummy);
	}
	double moving_peak::peak_function_hilly(const std::vector<double>& gen, int peak_number) {
		int j = 0;
		double dummy = (gen[0] - m_peak[peak_number][0])*(gen[0] - m_peak[peak_number][0]);
		for (j = 1; j < m_variable_size; j++)
			dummy += (gen[j] - m_peak[peak_number][j])*(gen[j] - m_peak[peak_number][j]);

		return m_height[peak_number] - m_width[peak_number] * dummy - 0.01*sin(20.0*dummy);
	}

	double moving_peak::peak_function_twin(const double * gen, int peak_number) {
		int j;
		double maximum = -100000.0, dummy = pow(gen[0] - m_peak[peak_number][0], 2);
		for (j = 1; j < m_variable_size; j++)
			dummy += pow(gen[j] - m_peak[peak_number][j], 2);

		dummy = m_height[peak_number] - m_width[peak_number] * dummy;

		maximum = dummy;
		dummy = pow(gen[0] - (m_peak[peak_number][0] + twin_peak[0]), 2);
		for (j = 1; j < m_variable_size; j++)
			dummy += pow(gen[j] - (m_peak[peak_number][j] + twin_peak[0]), 2);

		dummy = m_height[peak_number] + twin_peak[m_variable_size + 1] - ((m_width[peak_number] + twin_peak[m_variable_size])*dummy);
		if (dummy > maximum)
			maximum = dummy;

		return maximum;
	}

	double moving_peak::peak_function_twin(const std::vector<double>& gen, int peak_number) {
		int j;
		double maximum = -100000.0, dummy = pow(gen[0] - m_peak[peak_number][0], 2);
		for (j = 1; j < m_variable_size; j++)
			dummy += pow(gen[j] - m_peak[peak_number][j], 2);

		dummy = m_height[peak_number] - m_width[peak_number] * dummy;

		maximum = dummy;
		dummy = pow(gen[0] - (m_peak[peak_number][0] + twin_peak[0]), 2);
		for (j = 1; j < m_variable_size; j++)
			dummy += pow(gen[j] - (m_peak[peak_number][j] + twin_peak[0]), 2);

		dummy = m_height[peak_number] + twin_peak[m_variable_size + 1] - ((m_width[peak_number] + twin_peak[m_variable_size])*dummy);
		if (dummy > maximum)
			maximum = dummy;

		return maximum;
	}

	//TODO: can I delete this function?
	double moving_peak::function_selection(const double * gen, const int & peak_number) {
		double dummy = 0;
		switch (m_f) {
		case 1: {
			dummy = constant_basis_func(gen);
			break;
		}
		case 2: {
			dummy = five_peak_basis_func(gen);
			break;
		}
		case 3: {
			dummy = peak_function1(gen, peak_number);
			break;
		}
		case 4: {
			dummy = peak_function_cone(gen, peak_number);
			break;
		}
		case 5: {
			dummy = peak_function_hilly(gen, peak_number);
			break;
		}
		case 6: {
			dummy = peak_function_twin(gen, peak_number);
			break;
		}
		}
		return dummy;
	}
	double moving_peak::function_selection(const std::vector<double>& gen, const int & peak_number) {
		double dummy = 0;
		switch (m_f) {
		case 1: {
			dummy = constant_basis_func(gen);
			break;
		}
		case 2: {
			dummy = five_peak_basis_func(gen);
			break;
		}
		case 3: {
			dummy = peak_function1(gen, peak_number);
			break;
		}
		case 4: {
			dummy = peak_function_cone(gen, peak_number);
			break;
		}
		case 5: {
			dummy = peak_function_hilly(gen, peak_number);
			break;
		}
		case 6: {
			dummy = peak_function_twin(gen, peak_number);
			break;
		}
		}
		return dummy;
	}

	double moving_peak::dummy_eval(const double * gen) {
		int i;
		double maximum = -100000.0, dummy;

		for (i = 0; i < m_num_peaks; i++) {
			dummy = function_selection(gen, i);
			if (dummy > maximum)    maximum = dummy;
		}

		if (m_use_basis_function) {
			dummy = function_selection(gen, -1);
			/* If value of basis function is higher return it */
			if (maximum < dummy)    maximum = dummy;
		}
		return(maximum);
	}

	//TODO: need a set_record_flag in class optima.
	void moving_peak::initialize() {
		int i = 0;
		m_objective_accuracy =0.1;
		m_variable_accuracy = 0.2;
		set_range(0, 100);
		m_opt_mode[0] = optimization_mode::Maximization;
		m_optima.set_record_flag(false);
		update_time_linkage();
		for (i = 0; i < m_num_peaks; i++)
			for (int j = 0; j < m_variable_size; j++) {
				m_peak[i][j] = 100.0*global::ms_global->m_uniform[caller::Problem]->next();
				m_prev_movement[i][j] = global::ms_global->m_uniform[caller::Problem]->next() - 0.5;
			}

		if (m_standard_height <= 0.0) {
			for (i = 0; i < m_num_peaks; i++) m_height[i] = (m_max_height - m_min_height)*global::ms_global->m_uniform[caller::Problem]->next() + m_min_height;
		}
		else {
			for (i = 0; i < m_num_peaks; i++) m_height[i] = m_standard_height;
		}

		if (m_standard_width <= 0.0) {
			for (i = 0; i < m_num_peaks; i++)
				m_width[i] = (m_max_width - m_min_width)*global::ms_global->m_uniform[caller::Problem]->next() + m_min_width;
		}
		else {
			for (i = 0; i < m_num_peaks; i++)
				m_width[i] = m_standard_width;
		}

		calculate_global_optima();
		/*for (i=0; i< m_num_peaks; i++) {
		mp_heightOrder[i]=i;
		m_found[i]=false;
		}
		vector<int> idx(m_num_peaks);
		gQuickSort(m_height,m_num_peaks,idx);
		copy(idx.begin(),idx.end(),mp_heightOrder);
		gAmendSortedOrder<double*>(m_height,mp_heightOrder,mp_amendedHeightOrder,m_num_peaks);*/
		for (i = 0; i < m_num_peaks; i++) m_is_tracked[i] = 0;

		for (i = 0; i < m_num_peaks; i++)
			m_pre_peak[i] = m_peak[i];
		m_pre_height = m_height;
		m_pre_width = m_width;
		//calculate_associate_radius();
		m_peak_qaulity = 0;
		add_tag(problem_tag::MMP);
	}

	void moving_peak::current_peak_calculate(const double * gen) {
		int i;
		double maximum = -100000.0, dummy;

		m_current_peak = 0;
		maximum = function_selection(gen, 0);
		for (i = 1; i < m_num_peaks; i++) {
			dummy = function_selection(gen, i);
			if (dummy > maximum) {
				maximum = dummy;
				m_current_peak = i;
			}
		}
	}

	void moving_peak::allocate_memory(const int var_size, const int peaks) {
		m_shift.resize(var_size);
		m_covered_peaks.resize(peaks);
		m_prev_movement.resize(peaks);

		for (int i = 0; i < peaks; i++) {
			m_prev_movement.resize(var_size);
	}


}

	void moving_peak::random_change() {
		int i = 0, j = 0;
		double sum, sum2, offset;

		for (i = 0; i < m_num_peaks; i++) 	m_pre_peak[i] = m_peak[i];
		m_pre_height = m_height, m_height;
		m_pre_width = m_width;

		for (i = 0; i < m_num_peaks; i++) {
			if (m_whether_change[i] == false) continue;
			/* shift peak locations */
			sum = 0.0;
			for (j = 0; j < m_variable_size; j++) {
				m_shift[j] = global::ms_global->m_uniform[caller::Problem]->next() - 0.5;
				sum += m_shift[j] * m_shift[j];
			}
			if (sum > 0.0)		sum = m_vlength / sqrt(sum);
			else  sum = 0.0;                        /* only in case of rounding errors */

			sum2 = 0.0;
			for (j = 0; j < m_variable_size; j++) {
				m_shift[j] = sum*(1.0 - m_lambda)*m_shift[j] + m_lambda*m_prev_movement[i][j];
				sum2 += m_shift[j] * m_shift[j];
			}
			if (sum2 > 0.0)sum2 = m_vlength / sqrt(sum2);
			else     sum2 = 0.0;                      /* only in case of rounding errors */

			for (j = 0; j < m_variable_size; j++) {
				m_shift[j] *= sum2;
				m_prev_movement[i][j] = m_shift[j];
				if (m_domain[j].limit.first >(m_peak[i][j] + m_prev_movement[i][j])) {
					m_peak[i][j] = 2.0*m_domain[j].limit.first - m_peak[i][j] - m_prev_movement[i][j];
					m_prev_movement[i][j] *= -1.0;
				}
				else if (m_domain[j].limit.second< (m_peak[i][j] + m_prev_movement[i][j])) {
					m_peak[i][j] = 2.0*m_domain[j].limit.second- m_peak[i][j] - m_prev_movement[i][j];
					m_prev_movement[i][j] *= -1.0;
				}
				else
					m_peak[i][j] += m_prev_movement[i][j];
			}

			/* change peak width */

			offset = global::ms_global->m_normal[caller::Problem]->next()*m_width_severity[i];
			if ((m_width[i] + offset) < m_min_width)		m_width[i] = 2.0*m_min_width - m_width[i] - offset;
			else if ((m_width[i] + offset) > m_max_width)	m_width[i] = 2.0*m_max_width - m_width[i] - offset;
			else	m_width[i] += offset;

			if (change_counter() > 1 && m_change_peak_ratio < 1.0&&is_global_optima(i)) continue;
			/* change peak height */

			offset = m_height_severity[i] * global::ms_global->m_normal[caller::Problem]->next();

			if ((m_height[i] + offset) < m_min_height)	m_height[i] = 2.0*m_min_height - m_height[i] - offset;
			else if ((m_height[i] + offset) > m_max_height)	m_height[i] = 2.0*m_max_height - m_height[i] - offset;
			else	m_height[i] += offset;
		}

		calculate_global_optima();
		update_number_of_changes();

	}

	void moving_peak::change_num_peaks() {
		// TODO:need same() in class continuous
		// TODO:the constraint_value() in class problem can't be pure virtual
		unique_ptr<moving_peak> mpb = new moving_peak(m_variable_size, m_num_peaks_temp, m_change_peak_ratio, m_flag_variable_change
			, m_flag_num_peaks_change, m_num_peaks_change_mode, m_noise_flag, m_time_linkage_flag);

		mpb->copy(this);
		mpb->calculate_global_optima();

		allocate_memory(m_variable_size, m_num_peaks_temp);
		dynamic_continuous::allocate_memory(m_variable_size, m_num_peaks_temp);
		m_num_peaks = m_num_peaks_temp;
		*this = *mpb;
		mpb = 0;
	}

	//TODO: const parameter
	void moving_peak::copy(problem * problem) {
		dynamic_continuous::copy(problem);

		moving_peak *mpb = dynamic_cast<moving_peak *>(problem);
		unsigned dim = m_variable_size_temp < mpb->variable_size() ? m_variable_size_temp : mpb->variable_size();
		int peaks = m_num_peaks < mpb->number_of_peak() ? m_num_peaks : mpb->number_of_peak();

		m_f = mpb->m_f;
		m_vlength = mpb->m_vlength;
		m_lambda = mpb->m_lambda;
		m_use_basis_function = mpb->m_use_basis_function;
		m_calculate_right_peak = mpb->m_calculate_right_peak;
		m_standard_height = mpb->m_standard_height;
		m_standard_width = mpb->m_standard_width;

		//gCopy(m_shift,mpb->m_shift,dim);
		m_shift = mpb->m_shift;
		//gCopy(m_covered_peaks,mpb->m_covered_peaks,peaks);
		m_covered_peaks = mpb->m_covered_peaks;

		for (int i = 0; i < peaks; i++) {
			//gCopy(m_prev_movement[i],mpb->m_prev_movement[i],dim);
			m_prev_movement[i] = mpb->m_prev_movement[i];
		}
	}

	void moving_peak::update_peak_qaulity() {
		/*relative height over the heightest; Note perform just when a peak is found */
		m_peak_qaulity = 0;
		double sum = 0;
		for (int i = 0; i < m_num_peaks; i++) {
			if (m_found[i]) m_peak_qaulity += m_amended_height_order[i] * 1. / m_amended_height_order[m_height_order[m_num_peaks - 1]];
			sum += m_amended_height_order[i] * 1. / m_amended_height_order[m_height_order[m_num_peaks - 1]];
		}
		if (sum > 0)	m_peak_qaulity /= sum;
		else m_peak_qaulity = 0;
	}

	// TODO: should konw how to use the new evaluate()
	void moving_peak::calculate_associate_radius() {
		// to calculate an associate radius of peak i, find the nearest peak j, get the valley point, the distance
		// between peak i and the valley point is the associate radius of peak i
		vector<double> point(m_variable_size);
		vector<double> as_r(m_variable_size);

		for (int i = 0; i < m_num_peaks; i++) {
			double dis; int nearest = -1;
			if (!is_visable(i)) continue;
			for (int j = 0, count = 0; j < m_num_peaks; j++, count++) {
				if (j == i || !is_visable(j)) { count--; continue; }
				double d = 0;
				for (int dim = 0; dim < m_variable_size; dim++) {
					d += (m_peak[i][dim] - m_peak[j][dim])*(m_peak[i][dim] - m_peak[j][dim]);
				}
				d = sqrt(d);
				if (0 == count) {
					dis = d;
					nearest = j;
				}
				else if (dis > d) {
					dis = d;
					nearest = j;
				}
			}
			if (nearest != -1) {
				//normalize vector point
				for (int dim = 0; dim < m_variable_size; dim++) {
					point[dim] = m_peak[nearest][dim] - m_peak[i][dim];
					point[dim] /= dis;
				}


				double height, asHeight;
				height = asHeight = m_height[i];
				as_r = m_peak[i];
				//test in direction of point with a step of dis/100
				while (asHeight <= height) {
					bool flagBreak = false;
					for (int dim = 0; dim < m_variable_size; dim++) {
						as_r[dim] += dis / 100 * point[dim];
						if ((as_r[dim] - m_peak[i][dim])*(m_peak[nearest][dim] - as_r[dim]) < 0) {
							flagBreak = true;
							break;
						}
					}
					if (flagBreak) break;
					height = asHeight;
					solution<variable<real>,real> s(m_objective_size, m_variable_size);
					std::copy(as_r.begin(), as_r.end(), s.get_variable().begin());
					// TODO: should konw how to use the new evaluate()
					evaluate_(s, caller::Problem, false, true);
					asHeight = s.get_objective()[0];
				}

				m_associate_radius[i] = 0;
				for (int dim = 0; dim < m_variable_size; dim++) {
					// correction for one step backward
					as_r[dim] -= dis / 200 * point[dim];
					m_associate_radius[i] += (as_r[dim] - m_peak[i][dim])*(as_r[dim] - m_peak[i][dim]);
				}
				m_associate_radius[i] = sqrt(m_associate_radius[i]);
			}
			else {
				vector<double> r(2 * m_variable_size);
				for (int dim = 0; dim < m_variable_size; dim++) {
					double u, l;
					//HACK: no getSearchRange in class domain
					l = m_domain[dim].limit.first;
					u = m_domain[dim].limit.second;
					r[dim * 2] = fabs(l - m_peak[i][dim]);
					r[dim * 2 + 1] = fabs(u - m_peak[i][dim]);
				}
				m_associate_radius[i] = *min_element(r.begin(), r.end());
			}
		}
	}

	void moving_peak::set_severity() {
		for (int i = 0; i < m_num_peaks; i++) {
			auto get_rand_float = [](double min, double max)->double {return min + (max - min)*global::ms_global->m_uniform[caller::Problem]->next(); };
			m_height_severity[i] = get_rand_float(1, 10);//1.+9.*i/m_num_peaks;//severity of height changes, larger numbers  mean larger severity. in the contex of ROOT, peaks have different values
			m_width_severity[i] = get_rand_float(0.1, 1);// 0.1+0.9*i/m_num_peaks;severity of width changes, larger numbers mean larger severity
		}
	}

	moving_peak::moving_peak(const int var_size, const int num_peaks, double const changing_ratio, const bool flag_var_change, const bool flag_num_peak_change, const int peak_num_change_mode, const bool flag_noise, const bool flag_time_linkage) :\
		problem(string(), var_size, 1), dynamic_continuous(var_size, num_peaks, 1) {
		set_variable_change(flag_var_change);
		set_num_peak_change_mode(peak_num_change_mode);
		set_num_peaks_change(flag_num_peak_change);
		set_noise_flag(flag_noise);
		set_time_linkage_flag(flag_time_linkage);

		if (!readData()) exit(0);
		m_peaks_found = 0;

		allocate_memory(m_variable_size, m_num_peaks);
		m_name = "DYN_CONT_MovingPeak";
		initialize();
		set_number_of_changes((int)(changing_ratio*m_num_peaks));
		//TODO: setFlagLocTrue()?
		//m_optima.setFlagLocTrue();
		m_parameters << "Vlength:" << m_vlength << "; ";
	}

	moving_peak::moving_peak(param_map & v):problem(string(), (v[param_numDim]),  1), \
		dynamic_continuous((v[param_numDim]), (v[param_numPeak]), 1) {
		set_variable_change((v[param_flagNumDimChange]));
		set_num_peak_change_mode((v[param_peakNumChangeMode]));
		set_num_peaks_change((v[param_flagNumPeakChange]));
		set_noise_flag((v[param_flagNoise]));
		set_time_linkage_flag((v[param_flagTimeLinkage]));
		set_change_interval((v[param_changeFre]));
		set_change_type(change_type::CT_random);

		m_vlength = (v[param_shiftLength]);
		if ((v[param_flagNoise]) == true)		set_noise_severity_((v[param_noiseSeverity]));
		if ((v[param_flagTimeLinkage]) == true)	set_time_linkage_severity((v[param_timelinkageSeverity]));

		if (!readData()) exit(0);
		m_peaks_found = 0;
		allocate_memory(m_variable_size, m_num_peaks);
		m_name = "DYN_CONT_MovingPeak";
		initialize();
		set_number_of_changes((double)v[param_changeRatio]);
		//TODO: setFlagLocTrue()?
		//m_optima.setFlagLocTrue();
		m_parameters << "Vlength:" << m_vlength << "; ";
	}

	void moving_peak::change_stepsize_random() {
		m_vlength = global::ms_global->m_normal[caller::Problem]->next();
	}

	void moving_peak::change_stepsize_linear() {
		static	thread_local unique_ptr<int> counter;
		if (!counter.get()) counter.reset(new int(1));
		thread_local unique_ptr<double> frequency;
		if (!frequency.get()) frequency.reset(new double(3.14159 / 20.0));

		m_vlength = 1 + sin((double)(*counter)*(*frequency));
		(*counter)++;
	}

	int moving_peak::get_right_peak() {
		bool flag = false;
		for (int i = 0; i < m_num_peaks; i++) {
			if (m_optima_idx[i] == true && m_current_peak == i) {
				flag = true;
				break;
			}
		}
		return flag;
	}

	void moving_peak::set_v_length(const double s) {
		m_vlength = s;

		size_t start, end = 0;
		start = m_parameters.str().find("Vlength:");
		for (auto i = start; i < m_parameters.str().size(); i++) {
			if (m_parameters.str()[i] == ';') {
				end = i;
				break;
			}
		}
		stringstream ss;
		ss << "Vlength:" << m_vlength << "; ";
		string result = m_parameters.str();
		result.replace(start, end - start + 1, ss.str());
		m_parameters.str(result);
	}

	void moving_peak::reset() {
		set_severity();

		int i = 0;
		for (i = 0; i < m_num_peaks; i++)
			for (int j = 0; j < m_variable_size; j++) {
				m_peak[i][j] = 100.0*global::ms_global->m_uniform[caller::Problem]->next();
				m_prev_movement[i][j] = global::ms_global->m_uniform[caller::Problem]->next() - 0.5;
			}

		if (m_standard_height <= 0.0) {
			for (i = 0; i < m_num_peaks; i++) m_height[i] = (m_max_height - m_min_height)*global::ms_global->m_uniform[caller::Problem]->next() + m_min_height;
		}
		else {
			for (i = 0; i < m_num_peaks; i++) m_height[i] = m_standard_height;
		}

		if (m_standard_width <= 0.0) {
			for (i = 0; i < m_num_peaks; i++)
				m_width[i] = (m_max_width - m_min_width)*global::ms_global->m_uniform[caller::Problem]->next() + m_min_width;
		}
		else {
			for (i = 0; i < m_num_peaks; i++)
				m_width[i] = m_standard_width;
		}

		calculate_global_optima();
		m_change_counter = 0;
		for (i = 0; i < m_num_peaks; i++) m_found[i] = false;



		for (i = 0; i < m_num_peaks; i++) 	m_pre_peak[i] = m_peak[i];
		m_pre_height = m_height;
		m_pre_width = m_width;

		//calculate_associate_radius();
		/*for (i=0; i< m_num_peaks; i++) {
		mp_heightOrder[i]=i;
		m_found[i]=false;
		}
		vector<int> idx(m_num_peaks);
		gQuickSort(m_height,m_num_peaks,idx);
		copy(idx.begin(),idx.end(),mp_heightOrder);


		gAmendSortedOrder<double*>(m_height,mp_heightOrder,mp_amendedHeightOrder,m_num_peaks);*/
		m_peak_qaulity = 0;
		m_peaks_found = 0;
	}

	moving_peak & moving_peak::operator=(moving_peak & other) {
		if (this == &other) return *this;

		if (m_variable_size != other.m_variable_size || m_num_peaks != other.m_num_peaks) throw myexcept("Moving Peak assignment@moving_peak::operator=");
		dynamic_continuous::operator=(other);

		m_f = other.m_f;
		m_vlength = other.m_vlength;
		m_lambda = other.m_lambda;
		m_use_basis_function = other.m_use_basis_function;
		m_calculate_right_peak = other.m_calculate_right_peak;
		m_standard_height = other.m_standard_height;
		m_standard_width = other.m_standard_width;

		m_shift = other.m_shift;
		m_covered_peaks = other.m_covered_peaks;


		for (int i = 0; i < m_num_peaks; i++) {
			m_prev_movement[i] = other.m_prev_movement[i];

		}
		return *this;
	}

	double moving_peak::v_Length() {
		return m_vlength;
	}
	//TODO: complete it
	evaluation_tag moving_peak::evaluate_(base & ss, caller call, bool effective_fes, bool constructed) {
		solution<variable<real>,real> &s = dynamic_cast<solution<variable<real>, real> &>(ss);
		std:vector<real>x(m_variable_size);
		std::copy(s.get_variable().begin(), s.get_variable().end(), x.begin());
		if (this->m_noise_flag)	add_noise(x);

		double maximum = LONG_MIN, dummy;

		for (int i = 0; i < m_num_peaks; i++) {
			//if(maximum>mp_height[i]) continue; //optimization on the obj evaluation
			dummy = function_selection(x, i);
			if (dummy > maximum)      maximum = dummy;
		}

		if (m_use_basis_function) {
			dummy = function_selection(x, -1);
			/* If value of basis function is higher return it */
			if (maximum < dummy)     maximum = dummy;
		}
		s.get_objective[0] = maximum;

		//TODO: need a initialize_world_best() in class solution
		/*if (effective_fes&&m_effective_eval%m_change_interval == 0) {
			solution<variable<real>,real>::initilizeWB(s);
		}*/

		if (effective_fes&&is_tracked(x, s.get_objective())) update_peak_qaulity();
		if (effective_fes)    m_effective_eval++;
		bool flag;
#ifdef OFEC_CONSOLE
		if (global::ms_global->m_algorithm != nullptr)	flag = !global::ms_global->m_algorithm->terminating();
		else flag = true;
#endif
#ifdef OFEC_DEMON
		flag = true;
#endif
		if (effective_fes&&m_effective_eval%m_change_interval == 0 && flag) {

			//g_mutexStream.lock();
			//cout<<"The number of changes: "<<m_changeCounter<<endl;
			//g_mutexStream.unlock();
			//for(int i=0;i<m_num_peaks;i++)		printPeak(i);
			dynamic::change();
			//for(int i=0;i<m_num_peaks;i++)		printPeak(i);
			//getchar();
		}
		evaluation_tag rf = evaluation_tag::Normal;
		if (effective_fes) {
			if (global::ms_global->m_algorithm != nullptr) {
				if (global::ms_global->m_algorithm->terminating()) { rf = evaluation_tag::Terminate; }
				else if (global::ms_global->m_problem->is_tag(problem_tag::DOP)) {
					if (dynamic_cast<dynamic*>(global::ms_global->m_problem.get())->flag_time_linkage() && dynamic_cast<dynamic*>(global::ms_global->m_problem.get())->trigger_time_linkage()) {
						rf = evaluation_tag::Change_timelinkage;
					}
					if ((global::ms_global->m_problem->evaluations() + 1) % (dynamic_cast<dynamic*>(global::ms_global->m_problem.get())->change_interval()) == 0) {
						rf = evaluation_tag::Problem_change_next_eval;
					}
					if (global::ms_global->m_problem->evaluations() % (dynamic_cast<dynamic*>(global::ms_global->m_problem.get())->change_interval()) == 0) {
						if (dynamic_cast<dynamic*>(global::ms_global->m_problem.get())->flag_variable_change()) {
							rf = evaluation_tag::Change_dimension;
						}
						rf = evaluation_tag::Problem_change;
					}
				}
			}
		}
		return rf;
	}

