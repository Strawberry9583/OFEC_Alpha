#include "dynamic_continuous.h"
#include<fstream> 
#include<xutility> //TODO: the old version has an include.h
namespace OFEC {
	bool dynamic_continuous::is_global_opt(int idx) {
		return m_global_optima_idx[idx];
	}

	double dynamic_continuous::get_global_max() const {
		return m_global_optima;
	}

	void dynamic_continuous::print_fun(std::ofstream & out) {
		for (size_t i = 0; i < m_num_peaks; i++) {
			for (size_t j = 0; j < m_variable_size; j++)
				out << m_peak[i][j] << " "; // 传出内容似乎也并没有什么不妥
			out << endl;
		}
	}

	const vector<double>& dynamic_continuous::get_peak(const int p) const {
		if (p < 0 || p >= m_num_peaks) {
			throw myexcept("Please give right value of peak index [0,] @dynamic_continuous::get_peak");
		}
		return m_peak[p];
	}

	const vector<vector<double>>& dynamic_continuous::get_all_peaks() const {
		return m_peak;
	}

	double dynamic_continuous::get_peak_height(const int p) const {
		if (p < 0 || p >= m_num_peaks) {
			throw myexcept("Please give right value of peak index [0,] @dynamic_continuous::get_peak_height");
		}
		return m_height[p];
	}

	double dynamic_continuous::get_pre_peak_height(const int p) const {
		if (p < 0 || p >= m_num_peaks) {
			throw myexcept("Please give right value of peak index [0,] @dynamic_continuous::get_pre_peak_height");
		}
		return m_pre_height[p];
	}

	double dynamic_continuous::get_pre_peak_width(const int p) const {
		if (p < 0 || p >= m_num_peaks) {
			throw myexcept("Please give right value of peak index [0,] dynamic_continuous::get_pre_peak_width");
		}
		return m_pre_width[p];
	}

	const vector<double>& dynamic_continuous::get_height() const {
		return m_height;
	}

	const vector<double>& dynamic_continuous::get_pre_peak(const int p) const {
		if (p < 0 || p >= m_num_peaks) {
			throw myexcept("Please give right value of peak index [0,] @dynamic_continuous::get_pre_peak ");
		}
		return  m_pre_peak[p];
	}

	const vector<bool>& dynamic_continuous::get_global_optima_idx() const {
		return m_global_optima_idx;
	}

	int dynamic_continuous::get_number_of_global_opt_peak() const {
		return m_max_peaks_number;
	}

	void dynamic_continuous::set_number_of_changes(const int n) {
		if (n<1 || n>m_num_peaks) {
			throw myexcept("the number of changing peaks is invalid@dynamic_continuous::set_number_of_changes");
		}

		m_num_change_peaks = n;
		m_change_peak_ratio = n / (double)m_num_peaks;
		update_number_of_changes();

		size_t start, end;
		start = m_parameters.str().find("Changing peaks ratio:");
		for (auto i = start; i < m_parameters.str().size(); i++) {
			if (m_parameters.str()[i] == ';') {
				end = i;
				break;
			}
		}
		stringstream ss;
		ss << "Changing peaks ratio:" << m_change_peak_ratio << "; ";
		string result = m_parameters.str();
		result.replace(start, end - start + 1, ss.str());
		m_parameters.str(result);
	}

	void dynamic_continuous::set_number_of_changes(const double ratio) {
		if (ratio < 0 || ratio>1) {
			throw myexcept("the ratio of changing peaks is invalid@dynamic_continuous::set_number_of_changes");
		}

		m_change_peak_ratio = ratio;
		m_num_change_peaks = (int)(m_num_peaks*m_change_peak_ratio) > 1 ? (int)(m_num_peaks*m_change_peak_ratio) : 1;
		update_number_of_changes();

		size_t start, end;
		start = m_parameters.str().find("Changing peaks ratio:");
		for (auto i = start; i < m_parameters.str().size(); i++) {
			if (m_parameters.str()[i] == ';') {
				end = i;
				break;
			}
		}
		stringstream ss;
		ss << "Changing peaks ratio:" << m_change_peak_ratio << "; ";
		string result = m_parameters.str();
		result.replace(start, end - start + 1, ss.str());
		m_parameters.str(result);
	}

	void dynamic_continuous::set_height_severity(const double s) {
		for (size_t i = 0; i < m_num_peaks; i++) m_height_severity[i] = s;
	}

	void dynamic_continuous::set_width_severity(const double s) {
		for (size_t i = 0; i < m_num_peaks; i++) 	m_width_severity[i] = s;
	}

	dynamic_continuous & dynamic_continuous::operator=(const dynamic_continuous & dynamic_continuous_pro) {
		if (this == &dynamic_continuous_pro) return *this;

		if (m_variable_size != dynamic_continuous_pro.m_variable_size || m_num_peaks != dynamic_continuous_pro.m_num_peaks) return *this;

		dynamic::operator=(dynamic_continuous_pro);
		continuous::operator=(dynamic_continuous_pro);

		for (size_t i = 0; i < m_num_peaks; i++) {
			m_peak[i] = dynamic_continuous_pro.m_peak[i];
			m_pre_peak[i] = dynamic_continuous_pro.m_pre_peak[i];
			m_initial_peak[i] = dynamic_continuous_pro.m_initial_peak[i];
		}
		m_height = dynamic_continuous_pro.m_height;
		m_width = dynamic_continuous_pro.m_width;
		m_pre_height = dynamic_continuous_pro.m_pre_height;
		m_pre_width = dynamic_continuous_pro.m_pre_width;
		m_fit = dynamic_continuous_pro.m_fit;
		m_whether_change = dynamic_continuous_pro.m_whether_change;

		m_min_height = dynamic_continuous_pro.m_min_height;
		m_max_height = dynamic_continuous_pro.m_max_height;

		m_min_width = dynamic_continuous_pro.m_min_width;
		m_max_width = dynamic_continuous_pro.m_max_width;

		m_height_severity = dynamic_continuous_pro.m_height_severity;
		m_width_severity = dynamic_continuous_pro.m_width_severity;

		m_global_optima = dynamic_continuous_pro.m_global_optima;
		m_global_optima_idx = dynamic_continuous_pro.m_global_optima_idx;

		m_current_best = dynamic_continuous_pro.m_current_best;
		m_current_peak = dynamic_continuous_pro.m_current_peak;
		m_max_peaks_number = dynamic_continuous_pro.m_max_peaks_number;

		m_num_change_peaks = dynamic_continuous_pro.m_num_change_peaks;
		m_change_peak_ratio = dynamic_continuous_pro.m_change_peak_ratio;

		m_num_visable_peaks = dynamic_continuous_pro.m_num_visable_peaks;

		m_is_tracked, dynamic_continuous_pro.m_is_tracked;
		m_height_order, dynamic_continuous_pro.m_height_order;
		m_found, dynamic_continuous_pro.m_found;
		m_peaks_found = dynamic_continuous_pro.m_peaks_found;
		m_time_linkage, dynamic_continuous_pro.m_time_linkage;

		m_domain = dynamic_continuous_pro.m_domain;

		m_amended_height_order, dynamic_continuous_pro.m_amended_height_order;
		m_associate_radius, dynamic_continuous_pro.m_associate_radius;
		m_peak_qaulity = dynamic_continuous_pro.m_peak_qaulity;
		m_optima = dynamic_continuous_pro.m_optima;

		return *this;
	}

	void dynamic_continuous::set_height(const vector<double>& h) {
		m_height = h;
	}

	void dynamic_continuous::set_position(const vector<vector<double>>& p) {
		for (size_t i = 0; i < m_num_peaks; i++) {
			m_peak[i] = p[i];
			m_initial_peak[i] = p[i];
		}
	}

	void dynamic_continuous::set_width(const double w) {
		for (size_t i = 0; i < m_num_peaks; i++)
			m_width[i] = w;
	}

	void dynamic_continuous::print_peak(const int idx) {
		cout << "the " << idx << "th peak, height: " << m_height[idx] << " ass radius: " << m_associate_radius[idx] << " position:" << endl;

		for (int i = 0; i < m_variable_size; i++) {
			cout << m_peak[idx][i] << " ";
		}
		cout << endl;
	}

	void dynamic_continuous::print_peaks(ofstream & out) {
		for (int j = 0; j < m_num_peaks; j++) {
			for (int i = 0; i < m_variable_size; i++) {
				out << m_peak[j][i] << " ";
			}
			out << m_height[j] << endl;
		}

		/* for(int j=0;j<m_num_peaks;j++){
		out<<"set label \"p_{"<<j+1<<"}\" at ";
		for(int i=0;i<m_variable_size;i++){
		if(i+1<m_variable_size)
		out<<m_peak[j][i]+0.1<<", ";
		else out<<m_peak[j][i];
		}
		out<<endl;
		}*/
	}

	int dynamic_continuous::get_num_of_visable_peaks() {
		return m_num_visable_peaks;
	}

	bool dynamic_continuous::is_visable(const int idx) {

	}

	int dynamic_continuous::get_track_number(int idex) {
		return m_is_tracked[idex];
	}

	bool dynamic_continuous::is_tracked(vector<double>& gen, vector<double>& obj) {
		bool flag = false, movepeaks = false;
		for (int i = 0; i < m_num_peaks; i++) {
			double dis = 0, dis1 = fabs(obj[0] - m_height[i]);
			for (int j = 0; j < m_variable_size; j++) dis += (gen[j] - m_peak[i][j])*(gen[j] - m_peak[i][j]);
			dis = sqrt(dis);
			if (dis <= m_variable_accuracy&&dis1 <= m_objective_accuracy) {
				// peak[i] assumed to be found
				int j = 0;
				while (m_height_order[j++] != i&&j < m_num_peaks);
				if (!m_found[i]) {
					m_is_tracked[j - 1]++;
					m_found[i] = true;
					m_peaks_found++;
					flag = true;
				}

			}
			if (dis < m_time_linkage_severity) {
				// move peak[i] to a near random position when it was tracked
				if (m_time_linkage_flag) {
					move_peak(i);
					update_time_linkage();
					movepeaks = true;
					m_flag_trigger_time_linkage = true;
				}
			}
		}
		if (movepeaks) {
#ifdef DEMON_OFEC
			calculateSamplePoints();
#endif
		}
		return flag;
	}

	bool dynamic_continuous::is_tracked(double * gen, vector<double>& obj) {
		bool flag = false, movepeaks = false;
		for (int i = 0; i < m_num_peaks; i++) {
			double dis = 0, dis1 = fabs(obj[0] - m_height[i]);
			for (int j = 0; j < m_variable_size; j++) dis += (gen[j] - m_peak[i][j])*(gen[j] - m_peak[i][j]);
			dis = sqrt(dis);
			if (dis <= m_variable_accuracy&&dis1 <= m_objective_accuracy) {
				// peak[i] assumed to be found
				int j = 0;
				while (m_height_order[j++] != i&&j < m_num_peaks);
				if (!m_found[i]) {
					m_is_tracked[j - 1]++;
					m_found[i] = true;
					m_peaks_found++;
					flag = true;
				}
			}
			if (dis < m_time_linkage_severity) {
				// move peak[i] to a near random position when it was tracked
				if (m_time_linkage_flag) {
					move_peak(i);
					update_time_linkage();
					movepeaks = true;
					m_flag_trigger_time_linkage = true;
				}
			}
		}
		if (movepeaks) {
#ifdef DEMON_OFEC
			calculateSamplePoints();
#endif
		}
		return flag;
	}

	int dynamic_continuous::get_peaks_found() {
		return m_peaks_found;
	}

	double dynamic_continuous::get_associate_radius(int idx) {
		return m_associate_radius[idx];
	}

	double dynamic_continuous::get_peaks_traced_qaulity() {
		return m_peak_qaulity;
	}
	bool dynamic_continuous::is_global_optima_tracked() {
		// the global optimum is assumed to be tracked if any one of the global optima is tracked	
		for (int i = 0; i < m_num_peaks; i++) {
			if (m_global_optima_idx[i] && m_found[i]) return true;
		}
		return false;
	}

	// TODO: Need a class my_vector here
	const double * dynamic_continuous::get_nearest_peak(const vector<double>& p) {
		/*int nearest = 0;
		my_vector peak(global::ms_global->m_problem->variable_size(), m_peak[0]);
		double dis = peak.getDis(p);
		for (int i = 1; i<m_num_peaks; i++) {
			copy(m_peak[i], m_peak[i] + global::ms_global->m_problem->variable_size(), peak.begin());
			double d = peak.getDis(p);
			if (d<dis) {
				dis = d;
				nearest = i;
			}
		}
		return m_peak[nearest];*/
	}
	void dynamic_continuous::copy(problem * p) { //TODO: should I use a pointer or a reference?
		dynamic::copy(p);
		continuous::copy(p);

		dynamic_continuous *dcp = dynamic_cast<dynamic_continuous *>(p);

		int dim = m_variable_size_temp < p->variable_size() ? m_variable_size_temp : p->variable_size();
		int peaks = m_num_peaks < dcp->number_of_peak() ? m_num_peaks : dcp->number_of_peak();

		for (int i = 0; i < peaks; i++) {
			m_peak[i] = dcp->m_peak[i]; //TODO: Does it will cause some problem to use "=" rather than copy()?
			m_pre_peak[i] = dcp->m_pre_peak[i];
			m_initial_peak[i] = dcp->m_initial_peak[i];
		}
		m_height = dcp->m_height;
		m_width = dcp->m_width;
		m_pre_height = dcp->m_pre_height;
		m_pre_width = dcp->m_pre_width;
		m_fit = dcp->m_fit;
		m_whether_change = dcp->m_whether_change;

		m_min_height = dcp->m_min_height;
		m_max_height = dcp->m_max_height;

		m_min_width = dcp->m_min_width;
		m_max_width = dcp->m_max_width;

		m_height_severity = dcp->m_height_severity;
		m_width_severity = dcp->m_width_severity;

		m_global_optima = dcp->m_global_optima;
		m_global_optima_idx = dcp->m_global_optima_idx;

		m_current_best = dcp->m_current_best;
		m_current_peak = dcp->m_current_peak;
		m_max_peaks_number = dcp->m_max_peaks_number;

		m_change_peak_ratio = dcp->m_change_peak_ratio;
		m_num_change_peaks = (int)(m_change_peak_ratio*peaks) > 1 ? (int)(m_change_peak_ratio*peaks) : 1;//dcp->m_num_change_peaks;

		m_is_tracked = dcp->m_is_tracked;
		m_height_order = dcp->m_height_order;
		m_found= dcp->m_found;
		m_peaks_found = dcp->m_peaks_found;
		m_time_linkage= dcp->m_time_linkage;
		m_amended_height_order= dcp->m_amended_height_order;
		m_associate_radius= dcp->m_associate_radius;

		for (int j = 0; j < dim; j++) {
			m_domain[j] = dcp->m_domain[j];
		}

		m_peak_qaulity = dcp->m_peak_qaulity;
	}
	void dynamic_continuous::allocate_memory(const int size_var, const int num_peaks) {

	}
	// need a class my_vector here
	void dynamic_continuous::calculate_global_optima() {
		//if (m_opt_mode[0] == optimization_mode::Maximization) m_global_optima = *max_element(m_height.begin(), m_height.end());
		//else m_global_optima = *min_element(m_height.begin(), m_height.end());

		////TODO: need a clear function number in class optima
		////m_optima.clear();
		//m_max_peaks_number = 0;
		//double mindis = LONG_MAX;
		//for (int i = 0; i<m_num_peaks; i++) {
		//	m_global_optima_idx[i] = false;
		//	if (m_height[i] == m_global_optima) {
		//		for (int j = 0; j<m_num_peaks; ++j) {
		//			if (j == i) continue;
		//			my_vector s1(m_variable_size, m_peak[i]), s2(m_variable_size, m_peak[j]);

		//			double dis = s1.getDis(s2);
		//			if (mindis>dis) {
		//				mindis = dis;
		//			}
		//		}
		//		m_max_peaks_number++;
		//		m_global_optima_idx[i] = true;
		//		solution<variable<real>,real> s(m_variable_size, m_objective_size);
		//		std::copy(m_peak[i].begin(), m_peak[i].end(), s.get_variable().begin());
		//		s.get_objective()[0] = m_height[i];
		//		m_optima.append(s.get_variable());
		//	}
	}
	//TODO:need a shuffleIndex() in class global
	void dynamic_continuous::update_number_of_changes() {
		//if (m_num_change_peaks == m_num_peaks) {
		//	for (int i = 0; i<m_num_peaks; i++) m_whether_change[i] = true;
		//	return;
		//}
		//std::vector<int> a(m_num_peaks);

		//global::ms_global->shuffleIndex(a, m_num_peaks, caller::Problem);
		//// make sure the global optimum changes always
		//int gopt = 0;
		//for (int i = 0; i<m_num_peaks; i++) {
		//	if (m_global_optima_idx[i]) {
		//		gopt = i;
		//		break;
		//	}
		//}
		//int gidx;
		//for (int i = 0; i<m_num_peaks; i++) {
		//	if (a[i] == gopt) {
		//		gidx = i;
		//		break;
		//	}
		//}
		//int t = a[0];
		//a[0] = a[gidx];
		//a[gidx] = t;

		//for (int i = 0; i<m_num_peaks; i++) m_whether_change[i] = false;
		//for (int i = 0; i<m_num_change_peaks; i++) m_whether_change[a[i]] = true;

	}
	//TODO: the evaluate_() in class problem looks so different that I can't use
	void dynamic_continuous::compute_num_visable_peaks() {
		/*m_num_visable_peaks = m_num_peaks;
		for (int i = 0; i<m_num_peaks; i++) {
			solution<variable<real>,real> s(m_variable_size, m_objective_size);
			std::copy(m_peak[i].begin(), m_peak[i].end(), s.get_variable().begin());
			evaluate_(s,caller::Problem);
			double height = s.get_objective()[0];
			switch (m_opt_mode[0]) {
			case optimization_mode::Minimization:
				if (height<m_height[i]) m_num_visable_peaks--;
				break;
			case optimization_mode::Maximization:
				if (height>m_height[i]) m_num_visable_peaks--;
				break;
			}
		}*/
	}
	void dynamic_continuous::add_noise(double * x_) {
		for (int d = 0; d<m_variable_size; d++) {
			double x = x_[d];
			x += m_noise_severity_*global::ms_global->m_normal[caller::Problem]->next();
			if (m_domain[d].limit.second<x) x = m_domain[d].limit.second; //TODO: not intuitive enough
			if (m_domain[d].limit.first>x) x = m_domain[d].limit.first;
			x_[d] = x;
		}
	}
	void dynamic_continuous::update_time_linkage() {
		if (!m_time_linkage_flag) return;
		double range;
		for (int j = 0; j<m_variable_size; j++) {
			range = fabs(double(m_domain[j].limit.second - m_domain[j].limit.first));
			m_time_linkage[j] = 0.2*range*(global::ms_global->m_uniform[caller::Problem]->next() - 0.5);
		}
	}
	void dynamic_continuous::move_peak(const int idx) {
		if (idx<0 || idx >= m_num_peaks) throw myexcept("index out of boundary @ dynamic_continuous::move_peak(const int idx)");
		for (int d = 0; d<m_variable_size; d++) {
			double x = m_peak[idx][d];
			x += m_time_linkage[d];
			if (m_domain[d].limit.second<x) x = m_domain[d].limit.second;
			if (m_domain[d].limit.first>x)  x = m_domain[d].limit.first;
			m_peak[idx][d] = x;
		}
	}
}