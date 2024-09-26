#ifndef WRAPPER_FLOW_CUTTER_H
#define WRAPPER_FLOW_CUTTER_H

#include <string>
#include <vector>
#include <tuple>

void run_flow_cutter(int nb_nodes, const std::vector<int>& start_edges, const std::vector<int>& end_edges, 
	std::vector<std::vector<int>>& tree_bags,
	std::vector<int>& start_tree_edges,
	std::vector<int>& end_tree_edges,
	std::string method,
	int nb_runs = 1000, int random_seed = 0);

#endif //PACE_H