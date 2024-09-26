#ifndef TREE_DECOMPOSITION_H
#define TREE_DECOMPOSITION_H

#include "array_id_func.h"
#include "cell.h"
#include <string>
#include <ostream>

void store_tree_decompostion_of_order(ArrayIDIDFunc tail, ArrayIDIDFunc head, const ArrayIDIDFunc&order,
	std::vector<std::vector<int>>& nodes_in_bag, std::vector<int>& start_tree_edges, std::vector<int>& end_tree_edges);
void print_tree_decompostion_of_multilevel_partition(std::ostream&out, const ArrayIDIDFunc&tail, const ArrayIDIDFunc&head, const ArrayIDIDFunc&to_input_node_id, const std::vector<Cell>&cell_list);

#endif
