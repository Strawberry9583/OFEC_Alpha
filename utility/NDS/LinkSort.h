#ifndef LINKSORT_H
#define LINKSORT_H

#include <vector>

namespace NDS {
	void LinkSort(const std::vector<std::vector<double>>& data, std::vector<int>& rank, int& comp);
	struct LS_node
	{
		LS_node(const int value, LS_node* last = nullptr, LS_node* next = nullptr) : m_value(value), m_last(last), m_next(next) {}
		const int m_value;
		LS_node* m_last;
		LS_node* m_next;
	};
	class LS_list
	{
	public:
		LS_list() : m_begin(nullptr), m_end(nullptr) {}
		LS_node* push_back(const int value) {
			LS_node* new_node(new LS_node(value));
			if (m_begin == nullptr)
				m_begin = new_node;
			if (m_end == nullptr)
				m_end = new_node;
			else {
				new_node->m_last = m_end;
				m_end->m_next = new_node;
				m_end = new_node;
			}
			return new_node;
		}
		void erase(LS_node* node) {
			if (node == m_begin && node == m_end) {
				m_begin = nullptr;
				m_end = nullptr;
				delete node;
			}
			else if (node == m_begin) {
				node->m_next->m_last = nullptr;
				m_begin = node->m_next;
				delete node;
			}
			else if (node == m_end) {
				node->m_last->m_next = nullptr;
				m_end = node->m_last;
				delete node;
			}
			else {
				node->m_last->m_next = node->m_next;
				node->m_next->m_last = node->m_last;
				delete node;
			}
		}
		LS_node* begin() { return m_begin; }
		LS_node* end() { return m_end; }
	private:
		LS_node* m_begin;
		LS_node* m_end;
	};
}

#endif // !LINKSORT_H
