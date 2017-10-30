#include "real_DBG.h"

namespace OFEC {
	real_DBG::real_DBG(const int size_var, const int number_peaks, const int size_object) :dynamic_continuous(size_var, number_peaks, size_object) {
		//ctor 
		allocate_memory(m_variable_size, m_num_peaks);
		m_rotation_planes.reserve(0);

		m_max_height = 100;
		m_min_height = 10;
		m_max_width = 10;
		m_min_width = 1;
		set_noisy_severity(0.8);
		set_height_severity(5.0);
		set_width_severity(0.5);
		set_width(5); /// value in 1-10
	}
	bool real_DBG::set_period(const int period) {
		if (period < 1) return false;
		dynamic::set_period(period);

		m_rotation_planes.reserve(m_period);
		for (int i = 0; i < m_period; i++) {
			m_rotation_planes[i].reserve(m_num_peaks);
			for (int j = 0; j < m_num_peaks; j++)
				m_rotation_planes[i][j].reserve(m_variable_size);
		}
		return true;
	}
	real_DBG & real_DBG::operator=(const real_DBG & p) {
		if (this == &p) return *this;

		dynamic_continuous::operator=(p);
		m_rotation_matrix = p.m_rotation_matrix;

		m_prediction = p.m_prediction;

		if (period() != p.period()) {
			throw myexcept("The period must be the same@real_DBG::operator=");
		}

		for (int i = 0; i < m_period; i++) {
			for (int j = 0; j < m_num_peaks; j++)
				m_rotation_planes[i][j] = p.m_rotation_planes[i][j];
		}

		return *this;

	}
	void real_DBG::reset() {
		m_change_type.counter = 0;
		m_change_counter = 0;
		vector<double> t(m_num_peaks);

		for (int i = 0; i < m_num_peaks; i++) {
			if (m_change_type.type == change_type::CT_chaotic)
				t[i] = m_min_height + (m_max_height - m_min_height)*global::ms_global->m_uniform[caller::Problem]->next();
			else
				t[i] = 50;
		}
		set_height(t);
		t.reserve(0);


		vector<vector<double>> position(m_num_peaks);
		for (int i = 0; i < m_num_peaks; i++)
			position[i].reserve(m_variable_size);
		for (int i = 0; i < m_num_peaks; i++) {
			for (int j = 0; j < m_variable_size; j++) {
				position[i][j] = m_domain[j].limit.first + (m_domain[j].limit.second - m_domain[j].limit.first)*global::ms_global->m_uniform[caller::Problem]->next();
			}
		}
		set_position(position);


		for (int i = 0; i < m_num_peaks; i++) {
			position[i].clear();
			position[i].shrink_to_fit();
		}
		position.clear();
		position.shrink_to_fit();

		for (int i = 0; i < m_num_peaks; i++) m_pre_peak[i] = m_peak[i];
		m_pre_height = m_height;
		m_pre_width = m_width;

		calculate_global_optima();
	}
	void real_DBG::copy(const problem * problem) {
		dynamic_continuous::copy(problem);

		const real_DBG *r_dbg = dynamic_cast<const real_DBG *>(problem);
		int dim = m_variable_size_temp < problem->variable_size() ? m_variable_size_temp : problem->variable_size();
		int peaks = m_num_peaks < r_dbg->number_of_peak() ? m_num_peaks : r_dbg->number_of_peak();

		m_prediction = r_dbg->m_prediction;

		if (m_change_type.type == change_type::CT_recurrent || m_change_type.type == change_type::CT_recurrent_noisy) {
			for (int i = 0; i < r_dbg->m_period; i++) {
				if (m_change_type.counter <= i) break;
				for (int j = 0; j < peaks; j++) {
					if (dim == m_variable_size_temp) {// the number of dimensions decreases
						for (int m = 0, k = 0; k < dim; k++, m++)
							if (r_dbg->m_rotation_planes[i][j][m] == dim) { k--; continue; }
							else
								m_rotation_planes[i][j][k] = r_dbg->m_rotation_planes[i][j][m];

					}
					else
						std::copy(r_dbg->m_rotation_planes[i][j].begin(), r_dbg->m_rotation_planes[i][j].end(), m_rotation_planes[i][j]);
				}
			}
		}
	}
	double real_DBG::standard_change(const change_type T, const double min, const double max) {
		double step, sign;
		switch (T) {
		case change_type::CT_small_step:
			step = -1 + 2 * global::ms_global->m_uniform[caller::Problem]->next();
			step = m_alpha*step*(max - min);
			break;
		case change_type::CT_random:
			step = global::ms_global->m_normal[caller::Problem]->next();
			break;
		case change_type::CT_large_step:
			step = -1 + 2 * global::ms_global->m_uniform[caller::Problem]->next();
			if (step > 0)sign = 1;
			else if (step < 0) sign = -1;
			else sign = 0;
			step = (m_alpha*sign + (m_max_alpha - m_alpha)*step)*(max - min);
			break;
		case change_type::CT_recurrent:
		case change_type::CT_chaotic:
		case change_type::CT_recurrent_noisy:
			break;
		}
		return step;
	}
	void real_DBG::allocate_memory(const int size_var, const int peaks) {
		m_rotation_matrix.reserve(peaks);
		for (auto i = 0; i < peaks; i++) {
			m_rotation_matrix[i].resize(size_var, size_var);
		}
	}
	void real_DBG::initialize(const change_type type, const bool flag_var_change, const bool flag_number_peak_change, const int peak_number_change_mode, const bool flag_noise, const bool flag_time_linkage) {
		set_variable_change(flag_var_change);
		set_num_peak_change_mode(peak_number_change_mode);
		set_num_peaks_change(flag_number_peak_change);
		set_noise_flag(flag_noise);
		set_time_linkage_flag(flag_time_linkage);


		if (m_flag_variable_change) {
			m_change_type.type = change_type::CT_random;
			m_dir_variable_change = true;
		}
		else if (m_flag_num_peaks_change) {
			m_change_type.type = change_type::CT_random;
			m_dir_num_peaks_change = true;
		}
		else if (m_noise_flag || m_time_linkage_flag) {
			m_change_type.type = change_type::CT_random;

		}
		else {
			m_change_type.type = type;
		}
		m_change_type.counter = 0;

		if (m_change_type.type == change_type::CT_recurrent || m_change_type.type == change_type::CT_recurrent_noisy)      set_period(12);
		else      set_period(0);


		std::vector<double>t(m_num_peaks);

		set_choatic_constant(3.67);

		for (int i = 0; i < m_num_peaks; i++) {
			if (m_change_type.type == change_type::CT_chaotic)
				t[i] = m_min_height + (m_max_height - m_min_height)*global::ms_global->m_uniform[caller::Problem]->next();
			else
				t[i] = 50;
		}
		set_height(t);


		std::vector<std::vector<double>> position(m_num_peaks);
		for (int i = 0; i < m_num_peaks; i++)
			position[i].resize(m_variable_size);
		for (int i = 0; i < m_num_peaks; i++) {
			for (int j = 0; j < m_variable_size; j++) {
				position[i][j] = m_domain[j].limit.first + (m_domain[j].limit.second - m_domain[j].limit.first)*global::ms_global->m_uniform[caller::Problem]->next();
			}
		}
		set_position(position);

		for (int i = 0; i < m_num_peaks; i++) std::copy(m_peak[i].begin(), m_peak[i].end(), m_pre_peak[i]);
		std::copy(m_height.begin(), m_height.end(), m_pre_height);
		std::copy(m_width.begin(), m_width.end(), m_pre_width);
	}
	void real_DBG::restore_infor() {
		for (int i = 0; i < m_num_peaks; i++) {
			if (!m_whether_change[i]) {
				m_peak[i] = m_pre_peak[i];
				m_height[i] = m_pre_height[i];
				m_width[i] = m_pre_width[i];
			}
		}
	}
	void real_DBG::reinitialize() {
		if (m_num_peaks != *dynamic::ms_init_num_peaks) {
			m_num_peaks_temp = *dynamic::ms_init_num_peaks;
			change_num_peaks();
		}
		if (m_variable_size != *dynamic::ms_init_variable_size) {
			m_variable_size_temp = *dynamic::ms_init_variable_size;
			change_variable();
		}
	}
	void real_DBG::correct_solution(double * x) {
		for (int j = 0; j < m_variable_size; j++) {
			if (m_domain[j].limit.second < x[j])
				x[j] = m_domain[j].limit.second;
			else if (m_domain[j].limit.first > x[j])
				x[j] = m_domain[j].limit.first;
		}
	}
	void real_DBG::correct_solution(vector<double> x) {
		for (int j = 0; j < m_variable_size; j++) {
			if (m_domain[j].limit.second < x[j])
				x[j] = m_domain[j].limit.second;
			else if (m_domain[j].limit.first > x[j])
				x[j] = m_domain[j].limit.first;
		}
	}
	void real_DBG::height_standard_change() {
		double step;
		for (int i = 0; i < m_num_peaks; i++) {
			if (m_whether_change[i] == false) continue;
			step = m_height_severity[i] * standard_change(get_change_type(), m_min_height, m_max_height);
			m_height[i] = m_height[i] + step;
			if (m_height[i] > m_max_height || m_height[i] < m_min_height) m_height[i] = m_height[i] - step;

		}
	}
	void real_DBG::position_standard_change(double angle) {
		if (get_change_type() == change_type::CT_chaotic) {
			//HACK: original g_chaotic_value() in global is disappeared, so I moved it here.
			auto g_chaotic_value = [](const double x, const double min, const double max, const double chaotic_constant = 3.54) {
				if (min > max) return -1.0;
				double chaotic_value;
				chaotic_value = (x - min) / (max - min);
				chaotic_value = chaotic_constant*chaotic_value*(1 - chaotic_value);
				return (min + chaotic_value*(max - min));
			};
			for (int i = 0; i < m_num_peaks; i++) {
				if (m_whether_change[i] == false) continue;
				for (int j = 0; j < m_variable_size; j++)
					m_peak[i][j] = g_chaotic_value(m_peak[i][j], m_domain[j].limit.first, m_domain[j].limit.second, m_chaotic_constant);
			}
			return;
		}

		// for each basic function of dimension n(even number) , R=R(l1,l2)*R(l3,l4)*....*R(li-1,li), 0<=li<=n

		vector<int> d(m_variable_size);
		matrix I(m_variable_size, m_variable_size);
		vector<double> gene(m_variable_size);
		for (int i = 0; i < m_num_peaks; i++) {
			if ((get_change_type() == change_type::CT_recurrent || get_change_type() == change_type::CT_recurrent_noisy) && m_change_type.counter >= period()) {
				std::copy(m_rotation_planes[m_change_type.counter%period()][i].begin(), m_rotation_planes[m_change_type.counter%period()][i].end(), d);
			}
			else {
				// TODO: need a shuffleIndex in global or somewhere else
				shuffle_index(d, m_variable_size, global::ms_global->m_uniform[caller::Problem].get());
				if (get_change_type() == change_type::CT_recurrent || get_change_type() == change_type::CT_recurrent_noisy)
					std::copy(d.begin(), d.end(), m_rotation_planes[m_change_type.counter][i]);
			}

			if ((get_change_type() == change_type::CT_recurrent || get_change_type() == change_type::CT_recurrent_noisy) && m_change_type.counter%period() == 0)
				std::copy(m_initial_peak[i].begin(), m_initial_peak[i].end(), m_peak[i]);

			if (m_whether_change[i] == false) continue;

			I.identify();
			for (int j = 0; j + 1 < m_variable_size; j += 2) {
				if (get_change_type() == change_type::CT_small_step || get_change_type() == change_type::CT_large_step || get_change_type() == change_type::CT_random)
					angle = standard_change(get_change_type(), -OFEC_PI, OFEC_PI);
				//TODO: is this the original set_rotation_angle()?
				I.set_rotation_axes(d[j], d[j + 1], angle);
				if (j == 0) m_rotation_matrix[i] = I;
				else
					m_rotation_matrix[i] *= I;
			}
			matrix m(m_variable_size, 1);
			// TODO: the parameter number is not correct,k
			m.set_row(m_peak[i].data(),m_variable_size);
			m *= m_rotation_matrix[i];
			std::copy(m[0].begin(), m[0].end(), gene);
			correct_solution(gene);
			m_peak[i] = gene;
		}
	}
}
