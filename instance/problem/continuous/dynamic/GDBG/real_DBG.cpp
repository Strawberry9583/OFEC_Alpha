#include "real_DBG.h"

namespace OFEC {
	real_DBG::real_DBG(const int size_var, const int number_peaks, const int size_object):dynamic_continuous(size_var, number_peaks, size_object) {
		//ctor 
		allocateMemory(m_variable_size, m_num_peaks);
		mppp_rotationPlanes.reserve(0);

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
		if (period<1) return false;
		dynamic::set_period(period);

		mppp_rotationPlanes.reserve(m_period);
		for (int i = 0; i<m_period; i++) {
			mppp_rotationPlanes[i].reserve(m_num_peaks);
			for (int j = 0; j<m_num_peaks; j++)
				mppp_rotationPlanes[i][j].reserve(m_variable_size);
		}
		return true;
	}
	real_DBG & real_DBG::operator=(const real_DBG & p) {
		if (this == &p) return *this;

		dynamic_continuous::operator=(p);
		mp_rotationMatrix = p.mp_rotationMatrix;

		m_prediction = p.m_prediction;

		if (period() != p.period()) {
			throw myexcept("The period must be the same@real_DBG::operator=");
		}

		for (int i = 0; i<m_period; i++) {
			for (int j = 0; j<m_num_peaks; j++)
				mppp_rotationPlanes[i][j] = p.mppp_rotationPlanes[i][j];
		}

		return *this;

	}
	void real_DBG::reset() {
		m_change_type.counter = 0;
		m_change_counter = 0;
		vector<double> t(m_num_peaks);

		for (int i = 0; i<m_num_peaks; i++) {
			if (m_change_type.type == change_type:: CT_chaotic)
				t[i] = m_min_height + (m_max_height - m_min_height)*global::ms_global->m_uniform[caller::Problem]->next();
			else
				t[i] = 50;
		}
		set_height(t);
		t.reserve(0);


		vector<vector<double>> position(m_num_peaks);
		for (int i = 0; i<m_num_peaks; i++)
			position[i].reserve(m_variable_size);
		for (int i = 0; i<m_num_peaks; i++) {
			for (int j = 0; j<m_variable_size; j++) {
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

		for (int i = 0; i<m_num_peaks; i++) m_pre_peak[i] = m_peak[i];
		m_pre_height = m_height;
		m_pre_width = m_width;

		calculate_global_optima();
	}
	double real_DBG::standardChange(const change_type T, const double min, const double max) {
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
			if (step>0)sign = 1;
			else if (step<0) sign = -1;
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
	void real_DBG::allocateMemory(const int size_var, const int peaks) {
		mp_rotationMatrix.reserve(peaks);
		for (auto i = 0; i<peaks; i++) {
			mp_rotationMatrix[i].resize(size_var, size_var);
		}
	}
	void real_DBG::restoreInfor() {
		for (int i = 0; i<m_num_peaks; i++) {
			if (!m_whether_change[i]) {
				m_peak[i] = m_pre_peak[i];
				m_height[i] =m_pre_height[i];
				m_width[i] = m_pre_width[i];
			}
		}
	}
	void real_DBG::reinitialize() {
		if (m_num_peaks != *dynamic::ms_init_num_peaks) {
			m_num_peaks_temp = *dynamic::ms_init_num_peaks;
			changeNumPeaks();
		}
		if (m_variable_size != *dynamic::ms_init_variable_size) {
			m_variable_size_temp = *dynamic::ms_init_variable_size;
			changeDimension();
		}
	}
	void real_DBG::correctSolution(double * x) {
		for (int j = 0; j<m_variable_size; j++) {
			if (m_domain[j].limit.second<x[j])
				x[j] = m_domain[j].limit.second;
			else if (m_domain[j].limit.first>x[j])
				x[j] = m_domain[j].limit.first;
		}
	}
	void real_DBG::correctSolution(vector<double> x) {
		for (int j = 0; j<m_variable_size; j++) {
			if (m_domain[j].limit.second<x[j])
				x[j] = m_domain[j].limit.second;
			else if (m_domain[j].limit.first>x[j])
				x[j] = m_domain[j].limit.first;
		}
	}
	void real_DBG::heightStandardChange() {
		double step;
		for (int i = 0; i<m_num_peaks; i++) {
			if (m_whether_change[i] == false) continue;
			step = m_height_severity[i] * standardChange(get_change_type(), m_min_height, m_max_height);
			m_height[i] = m_height[i] + step;
			if (m_height[i]>m_max_height || m_height[i]<m_min_height) m_height[i] = m_height[i] - step;

		}
	}
	void real_DBG::positionStandardChange(double angle) {
		if (get_change_type() == change_type::CT_chaotic) {
			//HACK: original g_chaotic_value() in global is disappeared, so I moved it here.
			auto g_chaotic_value = [](const double x, const double min, const double max, const double chaotic_constant = 3.54) {
				if (min>max) return -1.0;
				double chaotic_value;
				chaotic_value = (x - min) / (max - min);
				chaotic_value = chaotic_constant*chaotic_value*(1 - chaotic_value);
				return (min + chaotic_value*(max - min));
			};
			for (int i = 0; i<m_num_peaks; i++) {
				if (m_whether_change[i] == false) continue;
				for (int j = 0; j < m_variable_size; j++)
					m_peak[i][j] = g_chaotic_value(m_peak[i][j], m_domain[j].limit.first, m_domain[j].limit.second, m_chaotic_constant);
			}
			return;
		}

		// for each basic function of dimension n(even number) , R=R(l1,l2)*R(l3,l4)*....*R(li-1,li), 0<=li<=n

		vector<int> d(m_variable_size);
		matrix I(m_variable_size, m_variable_size);
		double gene(m_variable_size);
		for (int i = 0; i<m_num_peaks; i++) {
			if ((get_change_type() == change_type::CT_recurrent || get_change_type() == change_type::CT_recurrent_noisy) && m_change_type.counter >= period()) {
				copy(mppp_rotationPlanes[m_change_type.counter%period()][i], mppp_rotationPlanes[m_change_type.counter%period()][i] + m_variable_size, d);
			}
			else {
				global::ms_global->shuffleIndex(d, m_variable_size, caller::Problem);
				if (get_change_type() == change_type::CT_recurrent || get_change_type() == change_type::CT_recurrent_noisy)
					copy(d, d + m_variable_size, mppp_rotationPlanes[m_change_type.counter][i]);
			}

			if ((get_change_type() == change_type::CT_recurrent || get_change_type() == change_type::CT_recurrent_noisy) && m_change_type.counter%period() == 0)
				copy(mpp_initialPeak[i], mpp_initialPeak[i] + m_variable_size, m_peak[i]);

			if (m_whether_change[i] == false) continue;

			I.identify();
			for (int j = 0; j + 1<m_variable_size; j += 2) {
				if (get_change_type() == change_type::CT_small_step || get_change_type() == change_type::CT_large_step || get_change_type() == change_type::CT_random)
					angle = standardChange(get_change_type(), -OFEC_PI, OFEC_PI);
				//TODO: is this the original set_rotation_angle()?
				I.set_rotation_axes(d[j], d[j + 1], angle);
				if (j == 0) mp_rotationMatrix[i] = I;
				else
					mp_rotationMatrix[i] *= I;
			}
			matrix m(m_variable_size, 1);
			// TODO: need a set_row() using vector<double> parameters.
			m.set_row(m_peak[i], m_variable_size);
			m *= mp_rotationMatrix[i];
			copy(m[0].begin(), m[0].end(), gene);
			correctSolution(gene);
			m_peak[i] = gene;
		}
	}
}
