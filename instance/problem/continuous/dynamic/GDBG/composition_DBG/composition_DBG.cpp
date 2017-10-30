#include "composition_DBG.h"
namespace OFEC {
	void composition_DBG::set_com_boundary() {
		for (int j = 0; j<m_variable_size; j++) {
			m_com_boundary[0][j].limit.first = 100.;
			m_com_boundary[0][j].limit.second = -100.;
			m_com_boundary[1][j].limit.first = 5.;
			m_com_boundary[1][j].limit.second = -5.;
			m_com_boundary[2][j].limit.first = 0.5;
			m_com_boundary[2][j].limit.second = -0.5;
			m_com_boundary[3][j].limit.first = 100.;
			m_com_boundary[3][j].limit.second = -100.;
			m_com_boundary[4][j].limit.first = 32.;
			m_com_boundary[4][j].limit.second = -32.;
		}
	}

	void composition_DBG::get_com_boundary(const function_tag &f, double &l, double &u, const int rDimIdx) const {
		switch (f) {
		case function_tag::Sphere:
			l = m_com_boundary[0][rDimIdx].limit.first;
			u = m_com_boundary[0][rDimIdx].limit.second;
			break;
		case function_tag::Rastrigin:
			l = m_com_boundary[1][rDimIdx].limit.first;
			u = m_com_boundary[1][rDimIdx].limit.second;
			break;
		case function_tag::Weierstrass:
			l = m_com_boundary[2][rDimIdx].limit.first;
			u = m_com_boundary[2][rDimIdx].limit.second;
			break;
		case function_tag::Griewank:
			l = m_com_boundary[3][rDimIdx].limit.first;
			u = m_com_boundary[3][rDimIdx].limit.second;
			break;
		case function_tag::Ackley:
			l = m_com_boundary[4][rDimIdx].limit.first;
			u = m_com_boundary[4][rDimIdx].limit.second;
			break;
		default:
			throw myexcept("No the function in the basic component fs@composition_DBG::get_com_boundary()");
			break;
		}
	}

	void composition_DBG::initialize(const change_type T, const composition_DBG_function_id rF, const double rChangingRatio, const bool flag_dim_change = false, \
		const bool flag_num_peak_change = false, const int peak_num_change_mode = 1, const bool flagNoise = false, const bool flag_timelinkage = false) {
		real_DBG::initialize(T, flag_dim_change, flag_num_peak_change, peak_num_change_mode, flagNoise, flag_timelinkage);
		set_opt_mode(optimization_mode::Minimization);
		set_number_of_changes(rChangingRatio);
		if (!ms_fun_id.get()) ms_fun_id.reset(new composition_DBG_function_id());
		*composition_DBG::ms_fun_id = rF;
		m_objective_accuracy = 1.0;
		m_variable_accuracy = 0.1;
		//TODO: need a set_flag_loc_true() in class optima
		//m_optima.setFlagLocTrue();
		std::vector<function_tag> basic_fun(m_num_peaks);

		switch (*composition_DBG::ms_fun_id) {
		case composition_DBG_function_id::COMDBG_SPHERE:
			for (int i = 0; i<m_num_peaks; i++) basic_fun[i] = function_tag::Sphere;
			break;
		case composition_DBG_function_id::COMDBG_RASTRIGIN:
			for (int i = 0; i<m_num_peaks; i++) basic_fun[i] = function_tag::Rastrigin;
			break;
		case composition_DBG_function_id::COMDBG_GRIEWANK:
			for (int i = 0; i<m_num_peaks; i++) basic_fun[i] = function_tag::Griewank;
			break;
		case composition_DBG_function_id::COMDBG_ACKLEY:
			for (int i = 0; i<m_num_peaks; i++) basic_fun[i] = function_tag::Ackley;
			break;
		case composition_DBG_function_id::COMDBG_HYBRID:
			basic_fun[0] = function_tag::Sphere;		basic_fun[5] = function_tag::Sphere;
			basic_fun[1] = function_tag::Rastrigin;		basic_fun[6] = function_tag::Rastrigin;
			basic_fun[2] = function_tag::Weierstrass;	basic_fun[7] = function_tag::Weierstrass;
			basic_fun[3] = function_tag::Griewank;		basic_fun[8] = function_tag::Griewank;
			basic_fun[4] = function_tag::Ackley;		basic_fun[9] = function_tag::Ackley;
			for (int i = 10; i<m_num_peaks; i++) basic_fun[i] = function_tag::Sphere;
			break;
		}
		set_basic_function(basic_fun.data());

		std::vector<double>t(m_num_peaks);
		for (int i = 0; i<m_num_peaks; i++)t[i] = 1.;
		set_coverge_sevrity(t.data());
		set_stretch_severity();
		set_rotation_matrix();

		matrix m(1, m_variable_size);
		std::vector<double> gene(m_variable_size);
		for (int i = 0; i<m_num_peaks; i++) {
			for (int j = 0; j<m_variable_size; j++) { // calculate the estimate max value of funciton i
				gene[j] = m_domain[j].limit.second;
				gene[j] /= m_stretch_severity[i];
			}
			m.set_row(gene.data(), m_variable_size);
			m *= m_rotation_matrix[i];
			std::copy(m[0].begin(), m[0].end(), gene.begin());
			correct_solution(m_com_function[i], gene.data());
			m_fmax[i] = select_fun(m_com_function[i], gene.data());
			if (m_fmax[i] == 0)   throw myexcept("the estimation max value must be greater not equal to 0@composition_DBG::initialize()");

		}

		calculate_global_optima();
		update_time_linkage();
	}
	composition_DBG::composition_DBG(const int rDimNumber, const int rNumPeaks, const change_type rT, const composition_DBG_function_id rF, \
		const double rChangingRatio, const bool flag_dim_change, const bool rFlagNumPeakChange, \
		const int peak_num_change_mode, const bool flagNoise, const bool flag_timelinkage) :problem(string(), rDimNumber, 1), real_DBG(rDimNumber, rNumPeaks), m_com_boundary(m_num_com_funs) {
		//TODO:delete it 
		allocate_memory(m_variable_size, m_num_peaks);
		m_name = "DYN_CONT_composition_DBG";
		set_com_boundary();
		m_height_normalize_severity = 2000.;

		initialize(rT, rF, rChangingRatio, flag_dim_change, rFlagNumPeakChange, peak_num_change_mode, flagNoise, flag_timelinkage);
	}

	composition_DBG::composition_DBG(param_map & v):problem(string(),v[param_numDim],1),\
		real_DBG(v[param_numDim], v[param_numPeak]), m_com_boundary(m_num_com_funs) {

		allocate_memory(m_variable_size, m_num_peaks);
		m_name = "DYN_CONT_composition_DBG";
		set_com_boundary();
		m_height_normalize_severity = 2000.;
		set_change_interval((v[param_changeFre]));
		initialize(static_cast<change_type>((int)(v[param_changeType])), static_cast<composition_DBG_function_id>((int)(v[param_comDBGFunID])), \
			(v[param_changeRatio]), (v[param_flagNumDimChange]), (v[param_flagNumPeakChange]), \
			(v[param_peakNumChangeMode]), (v[param_flagNoise]), (v[param_flagTimeLinkage]));
	}
	composition_DBG::~composition_DBG() {
		free_memory();
	}

	composition_DBG & composition_DBG::operator=(const composition_DBG &com) {
		if (this == &com) return *this;
		if (m_variable_size != com.m_variable_size || m_num_peaks != com.m_num_peaks) throw myexcept("composition_DBG assignment");
		real_DBG::operator =(com);

		set_com_boundary();
		std::copy(com.m_converge_severity.begin(), com.m_converge_severity.end(), m_converge_severity.begin());
		std::copy(com.m_stretch_severity.begin(), com.m_stretch_severity.end(), m_stretch_severity.begin());

		m_height_normalize_severity = com.m_height_normalize_severity;
		std::copy(com.m_com_function.begin(), com.m_com_function.end(), m_com_function.begin());
		std::copy(com.m_fmax.begin(), com.m_fmax.end(), m_fmax.begin());

		return *this;
	}

	void composition_DBG::set_rotation_matrix() {
		// for each basic function of dimension n(even number), R=R(l1,l2)*R(l3,l4)*....*R(ln-1,ln), 0<=li<=n
		matrix I(m_variable_size, m_variable_size);

		std::vector<int> d(m_variable_size);
		shuffle_index(d, m_variable_size, global::ms_global->m_uniform[caller::Problem].get());
		for (int i = 0; i<m_num_peaks; i++) {
			for (int j = 0; j + 1<m_variable_size; j += 2) {
				double angle = 2 * OFEC_PI*global::ms_global->m_uniform[caller::Problem]->next();				// random angle for rotation plane of d[j]-d[j+1] from d[j]th axis to d[j+1]th axis
				I.set_rotation_axes(d[j], d[j + 1], angle);
				if (j == 0)  m_rotation_matrix[i] = I;
				else	m_rotation_matrix[i] *= I;
				I.identify();
			}
		}
	}

	void composition_DBG::set_coverge_sevrity(const double * cs) {
		std::copy(cs, cs + m_num_peaks, m_converge_severity);
	}

	void composition_DBG::set_stretch_severity() {
		for (int i = 0; i<m_num_peaks; i++) {
			double l, u;
			get_com_boundary(m_com_function[i], l, u);
			m_stretch_severity[i] = m_converge_severity[i] * (m_domain[0].limit.second - m_domain[0].limit.first) / (u - l);
		}
	}

	void composition_DBG::set_basic_function(const function_tag * bf) {
		std::copy(bf, bf + m_num_peaks, m_com_function);
	}

	evaluation_tag composition_DBG::evaluate_(base & ss, caller call, bool effective_fes, bool constructed) {
		solution<variable<real>,real> &s = dynamic_cast<solution<variable<real>, real> &>(ss);
		vector<double> x(m_variable_size);
		std::copy(s.get_variable().begin(), s.get_variable().end(), x.begin());

		if (this->m_noise_flag)	add_noise(x);

		vector<double> width(m_num_peaks, 0), fit(m_num_peaks);
		for (int i = 0; i<m_num_peaks; i++) { // calculate weight for each function		
			for (int j = 0; j<m_variable_size; j++)
				width[i] += (x[j] - m_peak[i][j])*(x[j] - m_peak[i][j]);
			if (width[i] != 0)	width[i] = exp(-sqrt(width[i] / (2 * m_variable_size*m_converge_severity[i] * m_converge_severity[i])));
		}

		for (int i = 0; i<m_num_peaks; i++) { // calculate objective value for each function
			for (int j = 0; j<m_variable_size; j++)	// calculate the objective value of tranformation function i
				x[j] = (x[j] - m_peak[i][j]) / m_stretch_severity[i];//((1+fabs(m_peak[i][j]/mp_searchRange[j].m_upper))*
			matrix m(1, m_variable_size);
			m.set_row(x.data(), m_variable_size);
			m *= m_rotation_matrix[i];
			std::copy(m[0].begin(), m[0].end(), x.begin());
			correct_solution(m_com_function[i], x.data());
			fit[i] = select_fun(m_com_function[i], x.data());
			fit[i] = m_height_normalize_severity*fit[i] / fabs(m_fmax[i]);
			std::copy(s.get_variable().begin(), s.get_variable().end(), x.begin());
		}
		double sumw = 0, wmax;
		wmax = *max_element(width.begin(), width.end());
		for (int i = 0; i<m_num_peaks; i++)
			if (width[i] != wmax)
				width[i] = width[i] * (1 - pow(wmax, 10));
		for (int i = 0; i<m_num_peaks; i++)
			sumw += width[i];
		for (int i = 0; i<m_num_peaks; i++)
			width[i] /= sumw;
		double obj = 0;
		for (int i = 0; i<m_num_peaks; i++)
			obj += width[i] * (fit[i] + m_height[i]);
		s.get_objective()[0] = obj;
		// TODO: need a initializeWB(s) in class solution
		//TODO: maybe the static ms_minmax_objective in class problem
		if (effective_fes&&m_effective_eval%m_change_interval == 0) solution<variable<real>,real>::initilizeWB(s);

		if (effective_fes) {
			is_tracked(x.data(), s.get_objective());
			m_effective_eval++;
		}

		bool flag;
#ifdef OFEC_CONSOLE
		if (global::ms_global->m_algorithm != nullptr)	flag = !global::ms_global->m_algorithm->terminating();
		else flag = true;
#endif
#ifdef OFEC_DEMON
		flag = true;
#endif
		if (effective_fes&&m_effective_eval%m_change_interval == 0 && flag) {
			dynamic::change();
			if (m_time_linkage_flag) update_time_linkage();
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

	void composition_DBG::copy(const problem * rDP) {
		real_DBG::copy(rDP);
		const real_DBG *rP = dynamic_cast<const real_DBG*>(rDP);

		const composition_DBG* r_dbg = static_cast<const composition_DBG*>(rP);

		int peaks = m_num_peaks<rP->number_of_peak() ? m_num_peaks : rP->number_of_peak();

		std::copy(r_dbg->m_com_function.begin(), r_dbg->m_com_function.end(), m_com_function.begin());
		//	setRotationMatrix();
		std::copy(r_dbg->m_converge_severity.begin(), r_dbg->m_converge_severity.end(), m_converge_severity);
		std::copy(r_dbg->m_stretch_severity.begin(), r_dbg->m_stretch_severity.end(), m_stretch_severity);
	}

	void composition_DBG::free_memory() {
		m_stretch_severity.resize(0);
		m_stretch_severity.shrink_to_fit();
		
		m_converge_severity.resize(0);
		m_stretch_severity.shrink_to_fit();

		m_com_function.resize(0);
		m_com_function.shrink_to_fit();

		m_fmax.resize(0);
		m_com_function.shrink_to_fit();
	}

	void composition_DBG::random_change() {
		for (int i = 0; i<m_num_peaks; i++) std::copy(m_peak[i].begin(), m_peak[i].end(), m_pre_peak[i].begin());
		std::copy(m_height.begin(), m_height.end(), m_pre_height.begin());
		std::copy(m_width.begin(), m_width.end() + m_num_peaks, m_pre_width.begin());

		//change the global minimum value of each function
		height_standard_change();
		//change the position of global optimum of each function randomly
		position_standard_change(0);

		restore_infor();
		calculate_global_optima();
		update_number_of_changes();
		m_change_type.counter++;
	}

	void composition_DBG::small_step_change() {
		for (int i = 0; i<m_num_peaks; i++) std::copy(m_peak[i].begin(), m_peak[i].end(), m_pre_peak[i].begin());
		std::copy(m_height.begin(), m_height.end(), m_pre_height.begin());
		std::copy(m_width.begin(), m_width.end(), m_pre_width.begin());

		height_standard_change();
		position_standard_change(0);

		restore_infor();
		calculate_global_optima();
		update_number_of_changes();
		m_change_type.counter++;
	}

	void composition_DBG::large_step_change() {
		for (int i = 0; i<m_num_peaks; i++) std::copy(m_peak[i].begin(), m_peak[i].end(), m_pre_peak[i].begin());
		std::copy(m_height.begin(), m_height.end(), m_pre_height.begin());
		std::copy(m_width.begin(), m_width.end(), m_pre_width.begin());

		height_standard_change();
		position_standard_change(0);

		restore_infor();
		calculate_global_optima();
		update_number_of_changes();
		m_change_type.counter++;

	}

	void composition_DBG::recurrent_change() {
		for (int i = 0; i<m_num_peaks; i++) std::copy(m_peak[i].begin(), m_peak[i].end(), m_pre_peak[i].begin());
		std::copy(m_height.begin(), m_height.end(), m_pre_height.begin());
		std::copy(m_width.begin(), m_width.end(), m_pre_width.begin());
		double initial_angle;
		double height_range = m_max_height - m_min_height;

		for (int i = 0; i<m_num_peaks; i++) {
			if (m_whether_change[i] == false) continue;
			initial_angle = (double)m_period*i / m_num_peaks;
			m_height[i] = m_min_height + height_range*(sin(2 * OFEC_PI*(m_change_type.counter + initial_angle) / m_period) + 1) / 2.;
		}
		initial_angle = OFEC_PI*(sin(2 * OFEC_PI*m_change_type.counter / m_period) + 1) / 12.;
		position_standard_change(initial_angle);

		restore_infor();
		calculate_global_optima();
		update_number_of_changes();
		m_change_type.counter++;

	}

	void composition_DBG::chaotic_change() {
		for (int i = 0; i<m_num_peaks; i++) std::copy(m_peak[i].begin(), m_peak[i].end(), m_pre_peak[i].begin());
		std::copy(m_height.begin(), m_height.end(), m_pre_height.begin());
		std::copy(m_width.begin(), m_width.end(), m_pre_width.begin());

		for (int i = 0; i<m_num_peaks; i++) {
			if (m_whether_change[i] == false) continue;
			// need a g_chaotic_value in global, or I write to class real_DGB
			m_height[i] = chaotic_value(m_height[i], m_min_height, m_max_height);
		}

		position_standard_change(0);

		restore_infor();
		calculate_global_optima();
		update_number_of_changes();
		m_change_type.counter++;

	}

	void composition_DBG::recurrent_noisy_change() {
		for (int i = 0; i<m_num_peaks; i++) std::copy(m_peak[i].begin(), m_peak[i].end(), m_pre_peak[i].begin());
		std::copy(m_height.begin(), m_height.end(), m_pre_height.begin());
		std::copy(m_width.begin(), m_width.end(), m_pre_width.begin());

		double initial_angle;
		double height_range = m_max_height - m_min_height;

		double noisy;
		for (int i = 0; i<m_num_peaks; i++) {
			if (m_whether_change[i] == false) continue;
			initial_angle = (double)period()*i / m_num_peaks;

			m_height[i] = sin_value_noisy(m_change_type.counter, m_min_height, m_max_height, height_range, initial_angle, m_noisy_severity);
		}

		initial_angle = OFEC_PI*(sin(2 * OFEC_PI*(m_change_type.counter) / m_period) + 1) / 12.;
		noisy = m_noisy_severity*global::ms_global->m_normal[caller::Problem]->next();
		position_standard_change(initial_angle + noisy);

		restore_infor();
		calculate_global_optima();
		update_number_of_changes();
		m_change_type.counter++;

	}

	void composition_DBG::allocate_memory(const int var_size, const int rPeaks) {
		for (int i = 0; i<m_num_com_funs; i++) {
			m_com_boundary[i].resize(var_size); 
		}

		m_converge_severity.resize(rPeaks);
		m_stretch_severity.resize(rPeaks);
		m_com_function.resize(rPeaks);
		m_fmax.resize(rPeaks);
	}

	void composition_DBG::change_variable() {
		/// no need to preserve  previous information, e.g., positions, heights, width....

		composition_DBG* r_dbg = new composition_DBG(m_variable_size_temp, m_num_peaks, m_change_type.type, *ms_fun_id, m_change_peak_ratio, m_flag_variable_change,
			m_flag_num_peaks_change, m_num_peaks_change_mode, m_noise_flag, m_time_linkage_flag);


		if (m_change_type.type == change_type::CT_recurrent || m_change_type.type == change_type::CT_recurrent_noisy) {
			r_dbg->set_period(m_period);
		}

		r_dbg->copy(this);
		r_dbg->calculate_global_optima();

		free_memory();
		real_DBG::free_memory();
		dynamic_continuous::free_memory();
		//TODO: maybe no use
		// problem::free_memory();

		//TODO: check those allocate_memory()
		problem::allocate_memory(m_variable_size_temp);
		continuous::allocate_memory(m_variable_size_temp);
		dynamic_continuous::allocate_memory(m_variable_size_temp, m_num_peaks);
		real_DBG::allocate_memory(m_variable_size_temp, m_num_peaks);
		allocate_memory(m_variable_size_temp, m_num_peaks);

		m_variable_size = m_variable_size_temp;
		set_period(m_period);
		*this = *r_dbg;
		delete r_dbg;
		r_dbg = 0;
	}

	void composition_DBG::change_num_peaks() {
		composition_DBG* r_dbg = new composition_DBG(m_variable_size, m_num_peaks_temp, m_change_type.type, *ms_fun_id,
			m_change_peak_ratio, m_flag_variable_change, m_flag_num_peaks_change,
			m_num_peaks_change_mode, m_noise_flag, m_time_linkage_flag);

		if (m_change_type.type == change_type::CT_recurrent || m_change_type.type == change_type::CT_recurrent_noisy)
			r_dbg->set_period(m_period);
		r_dbg->copy(this);

		r_dbg->calculate_global_optima();

		free_memory();
		//TODO: check the free_memory()
		real_DBG::free_memory();
		dynamic_continuous::free_memory();

		//TODO:
		dynamic_continuous::allocate_memory(m_variable_size, m_num_peaks_temp);
		real_DBG::allocate_memory(m_variable_size, m_num_peaks_temp);
		allocate_memory(m_variable_size, m_num_peaks_temp);

		m_num_peaks = m_num_peaks_temp;
		set_period(m_period);
		*this = *r_dbg;
		delete r_dbg;
		r_dbg = 0;
	}

	void composition_DBG::correct_solution(const function_tag & f, double * x) {
		double l, u;
		get_com_boundary(f, l, u);
		for (int j = 0; j<m_variable_size; j++) {
			if (x[j]>u)  x[j] = u;
			else if (x[j]<l)  x[j] = l;
		}
	}

	double composition_DBG::f_sphere(double * x) {
		double fit = 0;
		for (int i = 0; i<m_variable_size; i++)
			fit += x[i] * x[i];
		return fit;
	}

	double composition_DBG::f_rastrigin(double * x) {
		double fit = 0;
		for (int i = 0; i<m_variable_size; i++)
			fit = fit + x[i] * x[i] - 10.*cos(2 * OFEC_PI*x[i]) + 10.;
		return fit;
	}

	double composition_DBG::f_weierstrass(double * x) {
		double a = 0.5, b = 3;
		int kmax = 20;
		double fit = 0, s = 0;
		for (int i = 0; i<m_variable_size; i++)
			for (int k = 0; k <= kmax; k++)
				fit += pow(a, k)*cos(2 * OFEC_PI*pow(b, k)*(x[i] + 0.5));
		for (int k = 0; k <= kmax; k++)
			s += pow(a, k)*cos(2 * OFEC_PI*pow(b, k)*0.5);
		s = s*m_variable_size;
		return fit - s;
	}

	double composition_DBG::f_griewank(double * x) {
		double s1 = 0, s2 = 1;
		for (int i = 0; i<m_variable_size; i++) {
			s1 += x[i] * x[i] / 4000.;
			s2 *= cos(x[i] / sqrt((double)(i + 1)));
		}
		return s1 - s2 + 1.;
	}

	double composition_DBG::f_ackley(double * x) {
		double fitness = 0;
		double s1 = 0, s2 = 0;
		for (int i = 0; i<m_variable_size; i++) {
			s1 += x[i] * x[i];
			s2 += cos(2 * OFEC_PI*x[i]);
		}
		fitness = -20 * exp(-0.2*sqrt(s1 / m_variable_size)) - exp(s2 / m_variable_size) + 20 + OFEC_E;
		return fitness;
	}

	double composition_DBG::select_fun(const function_tag & f, double * x) {
		double value;
		switch (f) {
		case function_tag::Sphere:
			value = f_sphere(x);
			break;
		case function_tag::Rastrigin:
			value = f_rastrigin(x);
			break;
		case function_tag::Weierstrass:
			value = f_weierstrass(x);
			break;
		case function_tag::Griewank:
			value = f_griewank(x);
			break;
		case function_tag::Ackley:
			value = f_ackley(x);
			break;
		default:
			break;
		}
		return value;
	}


}
