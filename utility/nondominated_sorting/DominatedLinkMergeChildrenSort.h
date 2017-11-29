#pragma once
#pragma once

#include"../random/newran.h"
#include<list>
#include<vector>
#include<queue>
#include<stack>
#include<iostream>
#include<algorithm>

#include"../nondominated_sorting/quick_sort.h"

//#include"DominanceGraphLinkSort.h"

using namespace std;
using namespace OFEC;
/*

Question(1): ID or pointer ?

*/


namespace NDS {





	class DLMC_sort {
		enum D_NodeState
		{
			DN_empty, DN_calculated, DN_deleted, DN_changed
		};
		struct L_Node {
			int m_id;
			L_Node * m_before;
			L_Node * m_after;
			bool m_isUsed;
			L_Node(int id = -1, bool used = false , L_Node* before = nullptr, L_Node * after = nullptr) :m_id(id),
				m_before(before), m_after(after), m_isUsed(used) {};
			
			
		};
		struct D_Node {
			int m_id;// the index of m_Dgraph
			int m_rank;
			vector<double> m_obj;
			vector<L_Node> m_objNode;
			vector<L_Node> m_mergeLinkNode;
			vector<bool> m_objVisited;
			//L_Node m_parents;
			//L_Node m_children;
			L_Node m_indexPosition;
		//	L_Node m_usedPosition;
			D_NodeState m_state;
			D_Node() :m_id(0), m_rank(0), m_state(DN_empty){};
			void init(int id, int objSize) {
			//	m_parents.m_after = m_children.m_after = end;
				m_indexPosition.m_isUsed = true;
				m_indexPosition.m_id = m_id = id;
				m_objNode.resize(objSize);
				for (auto & it : m_objNode) {
					it.m_id = id;
				}
				m_mergeLinkNode.resize(objSize);
				for (auto & it : m_mergeLinkNode) {
					it.m_id = id;
				}
				m_objVisited.resize(objSize);
			}
		};

	private:
		int m_nodeNum;
		int m_objNum;
		unsigned m_maxSize;
		bool m_isChange;
		const int m_firstRank = 0;
		vector<D_Node> m_Dgraph;
		vector<L_Node> m_parents;
		vector<L_Node> m_children;
		L_Node m_usedIndexHead;
		L_Node m_unUsedIndexHead;
		L_Node m_LNEnd;
		vector<vector<L_Node>> m_dominated;
		vector<L_Node> m_objectiveHead;
		int m_compareTimes;
		int m_offsetVisited = 10;
	public:
		//unsigned m_maxSize;// the maximum size of nodes to be deal with
		//the maximun size of nodes to be dealt with
		DLMC_sort(unsigned maxSize, unsigned objSize) :m_nodeNum(0), m_objNum(objSize), m_maxSize(maxSize),
			m_isChange(false), m_Dgraph(maxSize), m_parents(maxSize),m_children(maxSize), m_dominated(maxSize), m_objectiveHead(objSize), m_compareTimes(0)
		{
			initDataMember();
			initUnusedDgraphIndexLink(0, maxSize);
			m_offsetVisited = max(objSize + 5, unsigned(10));
		}


		void deleteData(vector<int>  delId) {

			sort(delId.begin(), delId.end(), isgreater<int, int>);
			int curId(m_nodeNum);
			auto itDelId(delId.begin());
			L_Node * pointer(&m_usedIndexHead);
			L_Node * childrenPointer(nullptr);
			queue<int> toDeletedId;
			while (pointer->m_after->m_isUsed) {
				pointer = pointer->m_after;
				if (--curId == *itDelId) {
					toDeletedId.push(pointer->m_id);
					m_Dgraph[pointer->m_id].m_state = DN_deleted;
					childrenPointer = &m_children[pointer->m_id];
					if (pointer->m_id == 88) {
						int stop = -1;
					}
					while (childrenPointer->m_after->m_isUsed) {
						childrenPointer = childrenPointer->m_after;
						deleteL_Node(&m_dominated[childrenPointer->m_id][pointer->m_id]);

					}
					if (++itDelId == delId.end())break;
				}

			}
			vector<pair<int, int>> sortedRank(m_nodeNum);
			pointer = &m_usedIndexHead;
			curId = m_nodeNum;
			while (pointer->m_after->m_isUsed) {
				pointer = pointer->m_after;
				sortedRank[--curId].first = m_Dgraph[pointer->m_id].m_rank;
				sortedRank[curId].second = pointer->m_id;
			}
			sort(sortedRank.begin(), sortedRank.end());
			vector<int> visited(m_maxSize, -1);
			queue<int> possibleSon;
			int visisingId(0);
			for (int i(m_nodeNum - 1); i >= 0; i--) {
				if (m_Dgraph[sortedRank[i].second].m_state != DN_deleted) {
					pointer = &m_children[sortedRank[i].second];
					while (pointer->m_after->m_isUsed) {
						pointer = pointer->m_after;

						if (m_Dgraph[pointer->m_id].m_state == DN_deleted) {
							possibleSon.push(pointer->m_id);
							deleteL_Node(&m_dominated[sortedRank[i].second][pointer->m_id]);
							deleteL_Node(&m_dominated[pointer->m_id][sortedRank[i].second]);

						}
						visited[pointer->m_id] = sortedRank[i].second;
					}
					while (!possibleSon.empty()) {
						visisingId = possibleSon.front();
						possibleSon.pop();
						pointer = &m_children[visisingId];
						while (pointer->m_after->m_isUsed) {
							pointer = pointer->m_after;
							if (visited[pointer->m_id] != sortedRank[i].second) {
								visited[pointer->m_id] = sortedRank[i].second;
								if (m_Dgraph[pointer->m_id].m_state == DN_deleted) {
									possibleSon.push(pointer->m_id);
								}
								else if (isDirectlyDominated(m_Dgraph[sortedRank[i].second].m_obj, &m_parents[pointer->m_id])) {
									insertL_Node(&m_children[sortedRank[i].second], &m_dominated[sortedRank[i].second][pointer->m_id]);
									insertL_Node(&m_parents[pointer->m_id], &m_dominated[pointer->m_id][sortedRank[i].second]);
								}
							}

						}
					}
				}
			}
			
			while (!toDeletedId.empty()) {
				visisingId = toDeletedId.front();
				toDeletedId.pop();
				deleteL_Node(&m_Dgraph[visisingId].m_indexPosition);
				insertL_Node(&m_unUsedIndexHead, &m_Dgraph[visisingId].m_indexPosition);
				m_Dgraph[visisingId].m_state = DN_empty;
				for (auto & it : m_Dgraph[visisingId].m_objNode) {
					deleteL_Node(&it);
				}
				clearLink(&m_parents[visisingId]);
				clearLink(&m_children[visisingId]);
			}

			m_nodeNum -= delId.size();
			m_isChange = true;
		}

		/*void deleteDataMO(vector<int> delId) {

			sort(delId.begin(), delId.end(), isgreater<int, int>);
			int curId(m_nodeNum);
			auto itDelId(delId.begin());
			L_Node * pointer(&m_usedIndexHead);
			L_Node * childrenPointer(nullptr);
			queue<int> toDeletedId;
			while (pointer->m_after->m_isUsed) {
				pointer = pointer->m_after;
				if (--curId == *itDelId) {
					toDeletedId.push(pointer->m_id);
					m_Dgraph[pointer->m_id].m_state = DN_deleted;
					childrenPointer = &m_children[pointer->m_id];
					if (pointer->m_id == 88) {
						int stop = -1;
					}
					while (childrenPointer->m_after->m_isUsed) {
						childrenPointer = childrenPointer->m_after;
						deleteL_Node(&m_dominated[childrenPointer->m_id][pointer->m_id]);

					}
					if (++itDelId == delId.end())break;
				}

			}
			vector<pair<int, int>> sortedRank(m_nodeNum);
			pointer = &m_usedIndexHead;
			curId = m_nodeNum;
			while (pointer->m_after->m_isUsed) {
				pointer = pointer->m_after;
				sortedRank[--curId].first = m_Dgraph[pointer->m_id].m_rank;
				sortedRank[curId].second = pointer->m_id;
			}
			sort(sortedRank.begin(), sortedRank.end());
			vector<int> visited(m_maxSize, -1);
			queue<int> possibleSon;
			int visisingId(0);
			int visitedContent(0);
			for (int i(m_nodeNum - 1); i >= 0; i--) {
				if (m_Dgraph[sortedRank[i].second].m_state != DN_deleted) {
					pointer = &m_children[sortedRank[i].second];
					while (pointer->m_after->m_isUsed) {
						pointer = pointer->m_after;

						if (m_Dgraph[pointer->m_id].m_state == DN_deleted) {
							possibleSon.push(pointer->m_id);
							deleteL_Node(&m_dominated[sortedRank[i].second][pointer->m_id]);
							deleteL_Node(&m_dominated[pointer->m_id][sortedRank[i].second]);

						}
						visited[pointer->m_id] = visitedContent;
					}
					while (!possibleSon.empty()) {
						visisingId = possibleSon.front();
						possibleSon.pop();
						pointer = &m_children[visisingId];
						while (pointer->m_after->m_isUsed) {
							pointer = pointer->m_after;
							if (visited[pointer->m_id] <visitedContent) {
								visited[pointer->m_id] = visitedContent;
								if (m_Dgraph[pointer->m_id].m_state == DN_deleted) {
									possibleSon.push(pointer->m_id);
								}
								else if (isDirectlyDominated(m_Dgraph[sortedRank[i].second].m_obj, &m_parents[pointer->m_id])) {
									insertL_Node(&m_children[sortedRank[i].second], &m_dominated[sortedRank[i].second][pointer->m_id]);
									insertL_Node(&m_parents[pointer->m_id], &m_dominated[pointer->m_id][sortedRank[i].second]);
								}
							}

						}
					}
				}
				visitedContent += m_offsetVisited;
			}

			while (!toDeletedId.empty()) {
				visisingId = toDeletedId.front();
				toDeletedId.pop();
				deleteL_Node(&m_Dgraph[visisingId].m_indexPosition);
				insertL_Node(&m_unUsedIndexHead, &m_Dgraph[visisingId].m_indexPosition);
				m_Dgraph[visisingId].m_state = DN_empty;
				for (auto & it : m_Dgraph[visisingId].m_objNode) {
					deleteL_Node(&it);
				}
				clearLink(&m_parents[visisingId]);
				clearLink(&m_children[visisingId]);
			}

			m_nodeNum -= delId.size();
			m_isChange = true;
		}*/
		void insertData(const vector<vector<double>> & data) {
			vector<L_Node> newObjectiveHead;
			initDataSort(data, newObjectiveHead);
			vector<L_Node> oldObjectiveHead(m_objNum);
			swap(m_objectiveHead, oldObjectiveHead);
			for (auto & it : m_objectiveHead) {
				it.m_after = &m_LNEnd;
			}
			vector<L_Node> newParents(m_maxSize);
			for (auto & it : newParents) {
				it.m_after = &m_LNEnd;
			}
			vector<L_Node> newChildren(m_maxSize);
			for (auto & it : newChildren) {
				it.m_after = &m_LNEnd;
			}
			L_Node * iter(&m_usedIndexHead);
			while (iter->m_after->m_isUsed) {
				iter = iter->m_after;
				for (int objIndex(0); objIndex < m_objNum; ++objIndex) {
					m_Dgraph[iter->m_id].m_objVisited[objIndex] = false;
				}
				m_Dgraph[iter->m_id].m_state = DN_changed;
			}
			vector<L_Node> newMergeLink(m_objNum);
			for (auto & it : newMergeLink) {
				it.m_after = &m_LNEnd;
			}			
			vector<L_Node> oldMergeLink(m_objNum);
			for (auto & it : oldMergeLink) {
				it.m_after = &m_LNEnd;
			}
			vector<L_Node*> oldObjectiveIter(m_objNum);
			for (int i(0); i < m_objNum; ++i) {
				oldObjectiveIter[i] = &oldObjectiveHead[i];
			}
			vector<L_Node*> newObjectiveIter(m_objNum);
			for (int i(0); i < m_objNum; ++i) {
				newObjectiveIter[i] = &newObjectiveHead[i];
			}
			vector<L_Node*> currentObjectiveIter(m_objNum);
			for (int i(0); i < m_objNum; ++i) {
				currentObjectiveIter[i] = &m_objectiveHead[i];
			}
			vector<int> visited(m_maxSize, -1);
			vector<int> compareIndex(m_maxSize, 0);
			vector<int> solvedId(m_nodeNum);
			int visitedContent(0);
			int solvedNode(0);
			L_Node * curPointer(nullptr);
			bool belongToOld(true);
			int nodeInOneLink(0);
			while (nodeInOneLink < m_nodeNum) {
				for (int objIndex(0); objIndex < m_objNum; ++objIndex) {
					if (oldObjectiveIter[objIndex]->m_after->m_isUsed&&newObjectiveIter[objIndex]->m_after->m_isUsed) {
						if (greater(m_Dgraph[oldObjectiveIter[objIndex]->m_after->m_id].m_obj, m_Dgraph[newObjectiveIter[objIndex]->m_after->m_id].m_obj, objIndex, m_compareTimes)) {
							curPointer = oldObjectiveIter[objIndex]->m_after;
							belongToOld = true;
						}
						else {
							curPointer = newObjectiveIter[objIndex]->m_after;
							belongToOld = false;
						}
					}
					else if (oldObjectiveIter[objIndex]->m_after->m_isUsed) {
						curPointer = oldObjectiveIter[objIndex]->m_after;
						belongToOld = true;
					}
					else{
						curPointer = newObjectiveIter[objIndex]->m_after;
						belongToOld = false;
					}
					if (m_Dgraph[curPointer->m_id].m_state != DN_calculated) {
						if (belongToOld) {
							findDirectChildren(&newMergeLink[objIndex], newParents, newChildren, curPointer->m_id, compareIndex, objIndex, visited, visitedContent);
						}
						else {
							findDirectChildren(&oldMergeLink[objIndex], newParents, newChildren, curPointer->m_id, compareIndex, objIndex, visited, visitedContent);
						}
						solvedId[m_nodeNum - solvedNode - 1] = curPointer->m_id;
							++solvedNode;
						visitedContent += m_offsetVisited;
						m_Dgraph[curPointer->m_id].m_state = DN_calculated;
					}
					deleteL_Node(curPointer);
					insertL_Node(currentObjectiveIter[objIndex], curPointer);
					currentObjectiveIter[objIndex] = currentObjectiveIter[objIndex]->m_after;
					if (solvedNode < m_nodeNum) {
						if (belongToOld) {
							mergeChildrenNode(&oldMergeLink[objIndex], curPointer->m_id, objIndex);
						}
						else {
							mergeChildrenNode(&newMergeLink[objIndex], curPointer->m_id, objIndex);
						}
						m_Dgraph[curPointer->m_id].m_objVisited[objIndex] = true;
					}

				}
				++nodeInOneLink;
			}
			for (auto & it : newMergeLink) {
				clearLink(&it);
			}
			for (auto & it : oldMergeLink) {
				clearLink(&it);
			}
			for (auto it : solvedId) {
				curPointer = &m_children[it];
				while (curPointer->m_after->m_isUsed) {
					curPointer = curPointer->m_after;
					if (!isDirectlyDominated(m_Dgraph[it].m_obj, &newParents[curPointer->m_id])) {
						deleteL_Node(curPointer);
						deleteL_Node(&m_dominated[curPointer->m_id][it]);
					}
				}
				while (newChildren[it].m_after->m_isUsed) {
					curPointer = newChildren[it].m_after;
					deleteL_Node(curPointer);
					if (isDirectlyDominated(m_Dgraph[it].m_obj, &newParents[curPointer->m_id])) {
						insertL_Node(&m_children[it], curPointer);
					}
					else {
						deleteL_Node(&m_dominated[curPointer->m_id][it]);
					}
				}
				while (newParents[it].m_after->m_isUsed) {
					curPointer = newParents[it].m_after;
					deleteL_Node(curPointer);
					insertL_Node(&m_parents[it], curPointer);
				}
			}
			m_isChange = true;
		}

		void initData(const vector<vector<double>> & data) {
			initDataSort(data, m_objectiveHead); 
			m_isChange = true;
		}

		void getRank(vector<int> & nodeRank) {
			if (m_isChange) {
				topoSort();
			}
			nodeRank.resize(m_nodeNum);
			int id(m_nodeNum);
			L_Node * pointer(&m_usedIndexHead);
			while (pointer->m_after->m_isUsed) {
				pointer = pointer->m_after;
				nodeRank[--id] = m_Dgraph[pointer->m_id].m_rank;
				if (id == 26) {
					int stop = -1;
				}
			}
		}

		
	private:

		void sortOjbectiveLink(int dataSize, const vector<vector<double>> & data, const vector<int>& allocatedId, vector<L_Node>& objectiveHead) {
			objectiveHead.resize(m_objNum);
			vector<int> index(dataSize);
			for (int objIndex(0); objIndex < m_objNum; ++objIndex) {
				objectiveHead[objIndex].m_after = &m_LNEnd;
				quick_sort(data, index, objIndex);
				for (int i(0); i < dataSize; ++i) {
					insertL_Node(&objectiveHead[objIndex], &m_Dgraph[allocatedId[index[i]]].m_objNode[objIndex]);
				}
			}
		}
		void clearLink(L_Node * head) {
			while (head->m_after->m_isUsed) {
				deleteL_Node(head->m_after);
			}
		}
		void findDirectChildren(L_Node * mergeLinkHead,vector<L_Node> & parents,vector<L_Node> & children,int curId,vector<int> & compareIndex,int objIndex,vector<int> & visited,int visitedContent) {
			    int ret(0);
				bool isFindSame(false);
				L_Node * dominatedPointer(nullptr);
				L_Node  non_dominatedPointer(-1, false, nullptr, &m_LNEnd);
				queue<int> queVisiting;
				while (mergeLinkHead->m_after->m_isUsed) {
					mergeLinkHead = mergeLinkHead->m_after;
					compareIndex[mergeLinkHead->m_id] = 0;
					ret = greaterEqual(m_Dgraph[curId].m_obj, m_Dgraph[mergeLinkHead->m_id].m_obj, compareIndex[mergeLinkHead->m_id], objIndex, m_compareTimes);
					if (ret == 1) {
						visited[mergeLinkHead->m_id] = visitedContent+1;
						insertL_Node(&children[curId], &m_dominated[curId][mergeLinkHead->m_id]);
						insertL_Node(&parents[mergeLinkHead->m_id], &m_dominated[mergeLinkHead->m_id][curId]);
					}
					else if (ret == 0) {
						isFindSame = true;
						dominatedPointer = &m_parents[mergeLinkHead->m_id];
						visited[mergeLinkHead->m_id] = visitedContent;
						while (dominatedPointer->m_after->m_isUsed) {
							dominatedPointer = dominatedPointer->m_after;
							insertL_Node(&parents[curId], &m_dominated[curId][dominatedPointer->m_id]);
							insertL_Node(&children[dominatedPointer->m_id], &m_dominated[dominatedPointer->m_id][curId]);
						}

						dominatedPointer = &m_children[mergeLinkHead->m_id];
						while (dominatedPointer->m_after->m_isUsed) {
							dominatedPointer = dominatedPointer->m_after;
							insertL_Node(&children[curId], &m_dominated[curId][dominatedPointer->m_id]);
							insertL_Node(&parents[dominatedPointer->m_id], &m_dominated[dominatedPointer->m_id][curId]);
						}
						break;
					}
					else {
						visited[mergeLinkHead->m_id] = visitedContent;
						queVisiting.push(mergeLinkHead->m_id);
					}
				}
				int visitingId(0);
				queue<int> possibleChildren;
				if (!isFindSame) {
					// visit all the non-dominated nodes
					while (!queVisiting.empty()) {
						visitingId = queVisiting.front();
						queVisiting.pop();
						dominatedPointer = &m_children[visitingId];
						while (dominatedPointer->m_after->m_isUsed) {
							dominatedPointer = dominatedPointer->m_after;
							if (visited[dominatedPointer->m_id] < visitedContent) {
								compareIndex[dominatedPointer->m_id] = compareIndex[visitingId];
								ret = greaterEqual(m_Dgraph[curId].m_obj, m_Dgraph[dominatedPointer->m_id].m_obj, compareIndex[dominatedPointer->m_id], objIndex, m_compareTimes);
								if (ret == 1) {
									visited[dominatedPointer->m_id] = visitedContent + 1;
									possibleChildren.push(dominatedPointer->m_id);
								}
								else {
									visited[dominatedPointer->m_id] = visitedContent;
									queVisiting.push(dominatedPointer->m_id);
								}
							}
							else {
								compareIndex[dominatedPointer->m_id] =max(compareIndex[visitingId], compareIndex[dominatedPointer->m_id]);
							}
						}
					}
					bool isDominated(true);
					// find the directly dominated nodes
					while (!possibleChildren.empty()) {
						visitingId = possibleChildren.front();
						possibleChildren.pop();
						isDominated = true;
						dominatedPointer = &m_parents[visitingId];
						while (dominatedPointer->m_after->m_isUsed) {
							dominatedPointer = dominatedPointer->m_after;
							if (m_Dgraph[dominatedPointer->m_id].m_objVisited[objIndex]&&visited[dominatedPointer->m_id] != visitedContent) {
								isDominated = false;
								break;
							}
						}
						if (isDominated) {
							insertL_Node(&children[curId], &m_dominated[curId][visitingId]);
							insertL_Node(&parents[visitingId], &m_dominated[visitingId][curId]);
						}
					}
				}

			
		}

		void initDataSort( const vector<vector<double>> & data,vector<L_Node> & objectiveHead) {
			int dataSize(data.size());
			vector<int> allocatedId(dataSize);
			insertDataMemory(data, dataSize, allocatedId);
			sortOjbectiveLink(dataSize, data, allocatedId, objectiveHead);
			for (auto & it : allocatedId) {
				for (int objIndex(0); objIndex < m_objNum; ++objIndex) {
					m_Dgraph[it].m_objVisited[objIndex] = false;
				}
			}

			vector<L_Node*> objectiveLinkPointer(m_objNum, nullptr);
			for (int objIndex(0); objIndex < m_objNum; ++objIndex) {
				objectiveLinkPointer[objIndex] = &objectiveHead[objIndex];
			}
			vector<L_Node> mergeLink(m_objNum);
			for (auto & it : mergeLink) {
				it.m_after = &m_LNEnd;
			}
			vector<int> visited(m_maxSize, -1);
			vector<int> compareIndex(m_maxSize, 0);
			int visitedContent(0);
			int solvedNode(0);
			while (solvedNode < dataSize) {
				for (int objIndex(0); objIndex < m_objNum; ++objIndex) {
					objectiveLinkPointer[objIndex] = objectiveLinkPointer[objIndex]->m_after;
					if (m_Dgraph[objectiveLinkPointer[objIndex]->m_id].m_state != DN_calculated) {
						findDirectChildren(&mergeLink[objIndex], m_parents, m_children, objectiveLinkPointer[objIndex]->m_id, compareIndex, objIndex, visited,visitedContent);
						if (++solvedNode == dataSize)
							break;
						visitedContent += m_offsetVisited;
						m_Dgraph[objectiveLinkPointer[objIndex]->m_id].m_state = DN_calculated;
					}
					mergeChildrenNode(&mergeLink[objIndex], objectiveLinkPointer[objIndex]->m_id, objIndex);
					m_Dgraph[objectiveLinkPointer[objIndex]->m_id].m_objVisited[objIndex] = true;
				}
			}
			for (auto & it : mergeLink) {
				clearLink(&it);
			}
			m_nodeNum += dataSize;
		}
		void mergeChildrenNode(L_Node * head, int curId,int objIndex) {
			insertL_Node(head, &m_Dgraph[curId].m_mergeLinkNode[objIndex]);
			head = head->m_after;
			while (head->m_after->m_isUsed) {
				head = head->m_after;
				if (m_dominated[curId][head->m_id].m_isUsed) {
					deleteL_Node(head);
				}
			}
		}
		int insertDataMemory(const vector<double> & data) {
			L_Node * curNode(m_unUsedIndexHead.m_after);
			deleteL_Node(curNode);
			insertL_Node(&m_usedIndexHead, curNode);
			m_Dgraph[curNode->m_id].m_obj = data;
			return curNode->m_id;
		}

		void insertDataMemory(const vector<vector<double>> & data, int dataSize,vector<int> & allocatedId) {
			allocatedId.resize(dataSize);
			for (int i(0); i < dataSize; ++i) {
				allocatedId[i] = insertDataMemory(data[i]);
			}
		}
		bool isDirectlyDominated(const vector<double> &ancester, L_Node * directParentHead) {
			int compareResult(0);
			while (directParentHead->m_after->m_isUsed) {
				directParentHead = directParentHead->m_after;
				compareResult = compareDominated(ancester, m_Dgraph[directParentHead->m_id].m_obj);
				if (compareResult == 1) {
					return false;
				}
				else if (compareResult == 0) {
					return true;
				}
				
			}
			return true;
		}
		bool isDirectlyDominated(const vector<double> &ancester, L_Node * directParentHead,int curId,vector<int> & visited,int visitedContent) {
			if (visited[curId] < visitedContent) {
				int compareResult(0);
				while (directParentHead->m_after->m_isUsed) {
					directParentHead = directParentHead->m_after;
					if (visited[directParentHead->m_id] < visitedContent) {
						compareResult = compareDominated(ancester, m_Dgraph[directParentHead->m_id].m_obj);
						if (compareResult == 1) {
							visited[directParentHead->m_id] = visitedContent;
							visited[curId] = visitedContent + 2;
							break;
						}
						else if (compareResult == 0) {
							visited[directParentHead->m_id] = visitedContent+1;
							visited[curId] = visitedContent;
							break;
						}
						else {
							visited[directParentHead->m_id] = visitedContent + 2;
						}
					}
					visited[curId] = visitedContent;
				}
			}
			if(visited[curId]==visitedContent)
			return true;
			else return false;
		}

		// time O(n+m)
		void topoSort() {
			vector<int> indegree(m_maxSize, 0);
			queue<int> queNode; 
			L_Node * iterPointer(&m_usedIndexHead);
			L_Node * dominatedPointer(nullptr);
			while (iterPointer->m_after->m_isUsed) {
				iterPointer = iterPointer->m_after;
				dominatedPointer = &m_parents[iterPointer->m_id];
				while (dominatedPointer->m_after->m_isUsed) {
					dominatedPointer = dominatedPointer->m_after;
					++indegree[iterPointer->m_id];
				}
				m_Dgraph[iterPointer->m_id].m_rank = m_firstRank;
				if (indegree[iterPointer->m_id] == 0) {
					queNode.push(iterPointer->m_id);
				}
			}
			int curId(0);
			while (!queNode.empty()) {
				curId = queNode.front();
				queNode.pop();
				iterPointer = &m_children[curId];
				while (iterPointer->m_after->m_isUsed) {
					iterPointer = iterPointer->m_after;
					if (--indegree[iterPointer->m_id] == 0) {
						m_Dgraph[iterPointer->m_id].m_rank = m_Dgraph[curId].m_rank + 1;
						queNode.push(iterPointer->m_id);
					}
				}
			}
			m_isChange = false;
		}
		// initialize

		
		void initUnusedDgraphIndexLink(int iter, int maxIter) {
			m_unUsedIndexHead.m_after = &m_Dgraph[iter].m_indexPosition;
			m_Dgraph[iter].m_indexPosition.m_before = &m_unUsedIndexHead;
			m_Dgraph[iter].m_indexPosition.m_after = &m_Dgraph[iter + 1].m_indexPosition;
			++iter;
			while (iter + 1 < maxIter) {
				m_Dgraph[iter].m_indexPosition.m_before = &m_Dgraph[iter - 1].m_indexPosition;
				m_Dgraph[iter].m_indexPosition.m_after = &m_Dgraph[iter + 1].m_indexPosition;
				++iter;
			}
			m_Dgraph[iter].m_indexPosition.m_before = &m_Dgraph[iter - 1].m_indexPosition;
			m_Dgraph[iter].m_indexPosition.m_after = &m_LNEnd;

		}
		// O ( maxSize * maxSize) 
		void initDataMember() {
			m_usedIndexHead.m_after = &m_LNEnd;
			m_LNEnd.m_isUsed = false;
			for (int i(0); i < m_maxSize; ++i) {
				m_Dgraph[i].init(i, m_objNum);
			}

			for (int i(0); i < m_maxSize; ++i) {
				m_dominated[i].resize(m_maxSize);
				for (int j(0); j < m_maxSize; ++j) {
					m_dominated[i][j].m_id =j;
				}
			}
			for (auto & it : m_objectiveHead) {
				it.m_after = &m_LNEnd;
			}
			for (auto & it : m_parents) {
				it.m_after = &m_LNEnd;
			}
			for (auto & it : m_children) {
				it.m_after = &m_LNEnd;
			}
		}


		// basic function

		// L_Node
		void insertL_Node(L_Node * head, L_Node * cur) {
			cur->m_after = head->m_after;
			cur->m_before = head;
			cur->m_after->m_before = head->m_after = cur;
			cur->m_isUsed = true;
		}

		void deleteL_Node(L_Node * cur) {
			cur->m_after->m_before = cur->m_before;
			cur->m_before->m_after = cur->m_after;
			cur->m_isUsed = false;
		}



		int compareDominated(const vector<double> & father, const vector<double> & son) {
			int ret = 0;
			for (int objIndex(0); objIndex < m_objNum; ++objIndex) {
				if (father[objIndex] > son[objIndex]) {
					return -1;
				}
				else if (father[objIndex] < son[objIndex]) {
					ret = 1;
				}
			}
			return ret;
		}

	};



}