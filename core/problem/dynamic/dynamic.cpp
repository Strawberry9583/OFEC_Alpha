#include"dynamic.h"

namespace OFEC {

	thread_local unique_ptr<int> dynamic::ms_init_num_peaks, dynamic::ms_init_num_variable, dynamic::ms_num_instance;

#ifdef OFEC_DEMON
#include "../../../ui/Buffer/Scene.h"
	extern unique_ptr<Scene> msp_buffer;
#endif


	void dynamic::set_variable_change(const bool flag) {
		m_flag_variable_change = flag;

		size_t start, end;
		start = m_parameters.str().find("VariableChange:");
		for (size_t i = start; i < m_parameters.str().size(); i++) {
			if (m_parameters.str()[i] == ';') {
				end = i;
				break;
			}
		}
		stringstream ss;
		ss << "VariableChange:" << m_flag_variable_change << "; ";
		string result = m_parameters.str();
		result.replace(start, end - start + 1, ss.str());
		m_parameters.str(result);
	}

	void dynamic::set_change_dirction(const bool flag) {
		m_dir_variable_change = flag;
	}

	dynamic::dynamic(const int size_var, const int num_peaks, const unsigned size_obj) :problem(string(), size_var, size_obj), m_change_counter(0)
		, m_variable_number_temp(size_var), m_num_peaks(num_peaks), m_num_peaks_temp(num_peaks), m_noise_flag(false), m_time_linkage_flag(false), m_flag_trigger_time_linkage(false) {

		m_change_interval = 5000;
		m_change_type.type = change_type::CT_random;
		m_change_type.counter = 0;
		m_period = 0;
		m_flag_variable_change = false;
		m_dir_variable_change = true;
		m_synchronize = true;
		m_noisy_severity = 0.8;

		m_alpha = 0.04;
		m_max_alpha = 0.1;
		m_chaotic_constant = 3.67; //in [3.57,4]

		m_flag_num_peaks_change = false;
		m_flag_num_peaks_change = true;
		m_num_peaks_change_mode = 1;
		m_noise_severity_ = 0.01;
		m_time_linkage_severity = 0.1;

		m_parameters << "Change frequency:" << m_change_interval << "; " << "TotalEvals:" << global::ms_arg[param_maxEvals] << "; " << "Peaks:" << m_num_peaks << "; " << "NumPeaksChange:" << m_flag_num_peaks_change << "-" << m_num_peaks_change_mode << "; " <<
			"NoisyEnvioronments:" << m_noise_flag << "; NoiseSeverity:" << m_noise_severity_ << "; TimeLinkageEnvironments:" << m_time_linkage_flag << "; TimeLinkageSeverity:" << m_time_linkage_severity << "; VariableChange:" << m_flag_variable_change << "; ";

		if (!ms_num_instance.get()) ms_num_instance.reset(new int(0));
		if (!ms_init_num_peaks.get()) ms_init_num_peaks.reset(new int);
		if (!ms_init_num_variable.get())ms_init_num_variable.reset(new int);
		(*ms_num_instance)++;
		if (*ms_num_instance == 1) {
			*ms_init_num_peaks = m_num_peaks;
			*ms_init_num_variable = m_variable_size;
		}
		add_tag(problem_tag::DOP);
	}

	dynamic::~dynamic() {
		//dtor
	}

	dynamic & dynamic::operator=(const dynamic & dynamic) {
		if (this == &dynamic) return *this;

		if (m_variable_size != dynamic.m_variable_size) {
			throw myexcept("The number of variables must be same!@dynamic::operator=");
		}
		if (m_change_type.type != dynamic.m_change_type.type) {
			throw myexcept("The change type must be same!@dynamic::operator=");
		}
		if (m_num_peaks != dynamic.m_num_peaks) {
			throw myexcept("The number of peaks must be same!@dynamic::operator=");
		}
		problem::operator=(dynamic);


		m_change_type.counter = dynamic.m_change_type.counter;
		m_change_interval = dynamic.m_change_interval;
		m_period = dynamic.m_period;
		m_flag_variable_change = dynamic.m_flag_variable_change;
		m_dir_variable_change = dynamic.m_dir_variable_change;
		m_synchronize = dynamic.m_synchronize;
		m_noisy_severity = dynamic.m_noisy_severity;

		m_alpha = dynamic.m_alpha;
		m_max_alpha = dynamic.m_max_alpha;
		m_chaotic_constant = dynamic.m_chaotic_constant;

		m_flag_num_peaks_change = dynamic.m_flag_num_peaks_change;
		m_flag_num_peaks_change = dynamic.m_flag_num_peaks_change;
		m_num_peaks_change_mode = dynamic.m_num_peaks_change_mode;

		m_noise_flag = dynamic.m_noise_flag;
		m_time_linkage_flag = dynamic.m_time_linkage_flag;

		m_noise_severity_ = dynamic.m_noise_severity_;
		m_time_linkage_severity = dynamic.m_time_linkage_severity;
		m_flag_trigger_time_linkage = dynamic.m_flag_trigger_time_linkage;
		return *this;
	}

	void dynamic::set_num_peak_change_mode(const int mode) {
		m_num_peaks_change_mode = mode;
	}

	int dynamic::get_num_peak_change_mode() {
		return m_num_peaks_change_mode;
	}

	void dynamic::set_noise_flag(const bool flag) {
		m_noise_flag = flag;

		size_t start, end;
		start = m_parameters.str().find("NoisyEnvioronments:");
		for (size_t i = start; i < m_parameters.str().size(); i++) {
			if (m_parameters.str()[i] == ';') {
				end = i;
				break;
			}
		}

		stringstream ss;
		ss << "NoisyEnvioronments:" << m_noise_flag << "; ";
		string result = m_parameters.str();
		result.replace(start, end - start + 1, ss.str());
		m_parameters.str(result);
	}

	int dynamic::get_number_of_peak() const {
		return m_num_peaks;
	}

	void dynamic::set_time_linkage_flag(const bool flag) {
		m_time_linkage_flag = flag;
		size_t start, end;
		start = m_parameters.str().find("TimeLinkageEnvironments:");
		for (size_t i = start; i < m_parameters.str().size(); i++) {
			if (m_parameters.str()[i] == ';') {
				end = i;
				break;
			}
		}

		stringstream ss;
		ss << "TimeLinkageEnvironments:" << m_time_linkage_flag << "; ";
		string result = m_parameters.str();
		result.replace(start, end - start + 1, ss.str());
		m_parameters.str(result);
	}

	void dynamic::change() {
		m_change_counter++;
		switch (get_change_type()) {
		case change_type::CT_random:
			random_change();
			break;
		case change_type::CT_recurrent:
			recurrent_change();
			break;
		case change_type::CT_recurrent_noisy:
			recurrent_noisy_change();
			break;
		case change_type::CT_small_step:
			small_step_change();
			break;
		case change_type::CT_large_step:
			large_step_change();
			break;
		case change_type::CT_chaotic:
			chaotic_change();
			break;
		default:
			break;
		}

		if (m_flag_variable_change) {

			if (m_variable_size == msc_min_variable_number)
				m_dir_variable_change = true;
			if (m_variable_size == msc_max_variable_number)
				m_dir_variable_change = false;

			if (m_dir_variable_change == true) {
				m_variable_number_temp += 1;
			}
			else {
				m_variable_number_temp -= 1;
			}
			change_variable();
		}

		if (m_flag_num_peaks_change) {
			if (m_num_peaks_change_mode == 1 || m_num_peaks_change_mode == 2) {
				if ((unsigned int)m_num_peaks >= msc_max_num_peaks - 1) m_dir_num_peaks_change = false;
				if ((unsigned int)m_num_peaks <= msc_min_num_peaks + 1) m_dir_num_peaks_change = true;
				int step = 0;

				if (m_name == "DYN_CONT_CompositionDBG") step = 2;
				else if (m_name == "DYN_CONT_RotationDBG") step = 2;
				else if (m_name == "DYN_CONT_MovingPeak") step = 2;
				else step = 2;

				if (m_num_peaks_change_mode == 2) {
					// TODO: too long
					step = static_cast<int>(step / 2 + (5 * step / 2 - step / 2)*global::ms_global->m_uniform.at(caller::Problem)->next());
				}

				if (m_dir_num_peaks_change == true) {
					if (m_num_peaks + step <= msc_max_num_peaks)		m_num_peaks_temp = m_num_peaks + step;
					else m_num_peaks_temp = msc_max_num_peaks;
				}
				else {
					if (m_num_peaks - step >= msc_min_num_peaks)		m_num_peaks_temp = m_num_peaks - step;
					else m_num_peaks_temp = msc_min_num_peaks;
				}
			}
			else {
				//random change
				//TODO: too long
				m_num_peaks_temp = static_cast<int>(msc_min_num_peaks + (msc_max_num_peaks - msc_min_num_peaks)*global::ms_global->m_uniform.at(caller::Problem)->next());
			}
			change_num_peaks();
		}

		//TODO: need measure module and a get_optima function in class continuous, and the run_id in class global
		//if (mSingleObj::getSingleObj() != nullptr) {
		//	vector<double> gOpt;
		//	if (global::ms_global->m_problem->getObjglobalOpt(gOpt)) {
		//		mSingleObj::getSingleObj()->addGOpt(global::ms_global->m_runId, gOpt[0]);
		//	}
		//	else {
		//		std::cout << "err" << endl;
		//	}
		//}

#ifdef OFEC_DEMON
		msp_buffer->updateFitnessLandsacpe_();
#endif

	}

	double dynamic::sin_value_noisy(const int x, const double min, const double max, const double amplitude, const double angle, const double noisy_severity) {
		double y;
		double noisy, t;
		y = min + amplitude*(sin(2 * OFEC_PI*(x + angle) / m_period) + 1) / 2.;
		noisy = noisy_severity*global::ms_global->m_normal.at(caller::Problem)->next();
		t = y + noisy;
		if (t > min&&t < max) y = t;
		else y = t - noisy;
		return y;
	}

	double dynamic::chaotic_step(const double x, const double min, const double max, const double scale) {
		if (min > max) return -1;
		double chaotic_value;
		chaotic_value = (x - min) / (max - min);
		chaotic_value = m_chaotic_constant*chaotic_value*(1 - chaotic_value);
		//return fabs((min+chaotic_value*(max-min)-x)* global::scale);
		return chaotic_value*scale;
	}

	bool dynamic::predict_change(const int evals_more) {
		int fre = get_change_fre();
		int evals = evaluations() % fre;
		if (evals + evals_more >= fre) return true;
		else return false;
	}

	void dynamic::set_noise_severity_(double value) {
		m_noise_severity_ = value;
		size_t start, end;
		start = m_parameters.str().find("NoiseSeverity:");
		for (size_t i = start; i < m_parameters.str().size(); i++) {
			if (m_parameters.str()[i] == ';') {
				end = i;
				break;
			}
		}
		stringstream ss;
		ss << "NoiseSeverity:" << m_noise_severity_ << "; ";
		string result = m_parameters.str();
		result.replace(start, end - start + 1, ss.str());
		m_parameters.str(result);

	}





	void dynamic::set_change_interval(const int change_interval) {
		if (change_interval > 0) m_change_interval = change_interval;
		else {
			throw myexcept("Change frequncy must be greater than 0 @dynamic::set_change_interval");
		}

		size_t start, end;
		start = m_parameters.str().find("Change frequency:");
		for (size_t i = start; i < m_parameters.str().size(); i++) {
			if (m_parameters.str()[i] == ';') {
				end = i;
				break;
			}
		}
		stringstream ss;
		ss << "Change frequency:" << m_change_interval << "; ";
		string result = m_parameters.str();
		result.replace(start, end - start + 1, ss.str());
		m_parameters.str(result);
	}

	bool dynamic::set_period(const int period) {
		if (period >= 0) m_period = period;
		else {
			throw myexcept("period must be positive@ dynamic::set_period");
		}
		return true;
	}

	void dynamic::set_change_type(const s_change_type & change_type) {
		m_change_type = change_type;
	}

	void dynamic::set_change_type(const change_type type) {
		m_change_type.type = type;
	}

	void dynamic::set_num_peaks_change(const bool peaks_change) {
		m_flag_num_peaks_change = peaks_change;
		size_t start, end;
		start = m_parameters.str().find("NumPeaksChange:");
		for (size_t i = start; i < m_parameters.str().size(); i++) {
			if (m_parameters.str()[i] == ';') {
				end = i;
				break;
			}
		}
	}

	void dynamic::set_synchronize(const bool flag) {
		m_dir_variable_change = flag;
	}

	void dynamic::set_noisy_severity(const double severity) {
		m_noisy_severity = severity;
	}

	void dynamic::set_timelinkage_severity(double value) {
		m_time_linkage_severity = value;

		size_t start, end;
		start = m_parameters.str().find("TimeLinkageSeverity:");
		for (size_t i = start; i < m_parameters.str().size(); i++) {
			if (m_parameters.str()[i] == ';') {
				end = i;
				break;
			}
		}
		stringstream ss;
		ss << "TimeLinkageSeverity:" << m_time_linkage_severity << "; ";
		string result = m_parameters.str();
		result.replace(start, end - start + 1, ss.str());
		m_parameters.str(result);
	}

	int dynamic::get_initial_num_peaks() {
		return *ms_init_num_peaks;
	}

	void dynamic::copy(problem * dynamic_problem) {
		problem::copy(dynamic_problem);

		dynamic *d_p = dynamic_cast<dynamic*>(dynamic_problem);
		m_change_type = d_p->m_change_type;
		m_change_interval = d_p->m_change_interval;
		m_period = d_p->m_period;
		m_flag_variable_change = d_p->m_flag_variable_change;
		m_dir_variable_change = d_p->m_dir_variable_change;
		m_synchronize = d_p->m_synchronize;
		m_noisy_severity = d_p->m_noisy_severity;
		m_alpha = d_p->m_alpha;
		m_max_alpha = d_p->m_max_alpha;
		m_chaotic_constant = d_p->m_chaotic_constant;

		m_flag_num_peaks_change = d_p->m_flag_num_peaks_change;
		m_dir_num_peaks_change = d_p->m_dir_num_peaks_change;
		m_num_peaks_change_mode = d_p->m_num_peaks_change_mode;
		m_noise_flag = d_p->m_noise_flag;
		m_time_linkage_flag = d_p->m_time_linkage_flag;
		m_noise_severity_ = d_p->m_noise_severity_;
		m_time_linkage_severity = d_p->m_time_linkage_severity;
		m_flag_trigger_time_linkage = d_p->m_flag_trigger_time_linkage;
	}

}
