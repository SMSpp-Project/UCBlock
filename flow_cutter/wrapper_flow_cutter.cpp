#include "wrapper_flow_cutter.h"
#include "id_func.h"
#include "list_graph.h"
#include "multi_arc.h"
#include "sort_arc.h"
#include "chain.h"
#include "union_find.h"
#include "node_flow_cutter.h"
#include "separator.h"
#include "id_multi_func.h"
#include "filter.h"
#include "preorder.h"
#include "contraction_graph.h"
#include "tree_decomposition.h"
#include "greedy_order.h"
#include "min_max.h"

#include <limits>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sstream>
#include <queue>

#include <iostream>
using namespace std;

ArrayIDIDFunc tail, head;

// This hack is actually standard compilant
template <class T, class S, class C>
S& access_internal_vector(std::priority_queue<T, S, C>& q) {
	struct Hacked : private priority_queue<T, S, C> {
		static S& access(priority_queue<T, S, C>& q) {
			return q.*&Hacked::c;
		}
	};
	return Hacked::access(q);
}


template<class Tail, class Head>
void check_multilevel_partition_invariants(const Tail&tail, const Head&head, const std::vector<Cell>&multilevel_partition){
	#ifndef NDEBUG
	const int node_count = tail.image_count();
	const int arc_count = tail.preimage_count();

	auto is_child_of = [&](int c, int p){
		for(;;){
			if(c == p)
				return true;
			if(c == -1)
				return false;
			c = multilevel_partition[c].parent_cell;
		}
	};

	auto are_ordered = [&](int a, int b){
		return is_child_of(a, b) || is_child_of(b, a);
	};

	ArrayIDFunc<int> cell_of_node(node_count);
	cell_of_node.fill(-1);

	for(int i=0; i<(int)multilevel_partition.size(); ++i){
		for(auto&y:multilevel_partition[i].separator_node_list){
			assert(cell_of_node(y) == -1);
			cell_of_node[y] = i;
		}
	}

	for(auto x:cell_of_node)
		assert(x != -1);

	for(int xy = 0; xy < arc_count; ++xy){
		int x = cell_of_node(tail(xy)), y = cell_of_node(head(xy));
		assert(are_ordered(x, y));
	}
	#endif
}

template<class Tail, class Head, class ComputeSeparator, class OnNewMP>
void compute_multilevel_partition(const Tail&tail, const Head&head, const ComputeSeparator&compute_separator, int smallest_known_treewidth, const OnNewMP&on_new_multilevel_partition){

	const int node_count = tail.image_count();
	const int arc_count = tail.preimage_count();

	std::vector<Cell>closed_cells;
	std::priority_queue<Cell>open_cells;

	{
		Cell top_level_cell;
		top_level_cell.separator_node_list.resize(node_count);
		for(int i=0; i<node_count; ++i)
			top_level_cell.separator_node_list[i] = i;
		//top_level_cell.boundary_node_list = {};
		top_level_cell.parent_cell = -1;

		open_cells.push(std::move(top_level_cell));
	}

	int max_closed_bag_size = 0;
	int max_open_bag_size = node_count;

	auto check_if_better = [&]{
		int current_tree_width = std::max(max_closed_bag_size, max_open_bag_size);

		if(current_tree_width < smallest_known_treewidth){
			smallest_known_treewidth = current_tree_width;

			std::vector<Cell>cells = closed_cells;
			for(auto&q:access_internal_vector(open_cells))
				cells.push_back(q);
			check_multilevel_partition_invariants(tail, head, cells);
			on_new_multilevel_partition(cells, open_cells.empty() || max_closed_bag_size>=max_open_bag_size);
		}
	};

	check_if_better();


	ArrayIDFunc<int>node_to_sub_node(node_count);
	node_to_sub_node.fill(-1);

	auto inv_tail = invert_sorted_id_id_func(tail);

	BitIDFunc in_child_cell(node_count);
	in_child_cell.fill(false);

	while(!open_cells.empty()){

		#ifndef NDEBUG
		
		int real_max_closed_bag_size = 0;
		for(auto&x:closed_cells)
			max_to(real_max_closed_bag_size, x.bag_size());
		assert(max_closed_bag_size == real_max_closed_bag_size);

		int real_max_open_bag_size = 0;
		for(auto&x:access_internal_vector(open_cells))
			max_to(real_max_open_bag_size, x.bag_size());
		assert(max_open_bag_size == real_max_open_bag_size);

		#endif

		auto current_cell = std::move(open_cells.top());
		open_cells.pop();

		bool must_recompute_max_open_bag_size = (current_cell.bag_size() == max_open_bag_size);

		int closed_cell_id = closed_cells.size();

		if(current_cell.bag_size() > max_closed_bag_size){

			auto interior_node_list = std::move(current_cell.separator_node_list);
			int interior_node_count = interior_node_list.size();			

			ArrayIDFunc<int>sub_node_to_node(interior_node_count);

			int next_sub_id = 0;
			for(int x:interior_node_list){
				node_to_sub_node[x] = next_sub_id;
				sub_node_to_node[next_sub_id] = x;
				++next_sub_id;
			}

			auto is_node_interior = id_func(
				node_count,
				[&](int x)->bool{
					return node_to_sub_node(x) != -1;
				}
			);

			auto is_arc_interior = id_func(
				arc_count,
				[&](int xy)->bool{
					return is_node_interior(tail(xy)) && is_node_interior(head(xy));
				}
			);

			int interior_arc_count = count_true(is_arc_interior);
			auto sub_tail = keep_if(is_arc_interior, interior_arc_count, tail);
			auto sub_head = keep_if(is_arc_interior, interior_arc_count, head);

			for(auto&x:sub_tail)
				x = node_to_sub_node(x);
			sub_tail.set_image_count(interior_node_count);

			for(auto&x:sub_head)
				x = node_to_sub_node(x);
			sub_head.set_image_count(interior_node_count);
			
			auto sub_separator = compute_separator(sub_tail, sub_head);

			BitIDFunc is_in_sub_separator(interior_node_count);
			is_in_sub_separator.fill(false);
			for(auto x:sub_separator)
				is_in_sub_separator.set(x, true);

			UnionFind uf(interior_node_count);

			for(int xy=0; xy<interior_arc_count; ++xy){
				int x = sub_tail(xy);
				int y = sub_head(xy);
				if(!is_in_sub_separator(x) && !is_in_sub_separator(y))
					uf.unite(x, y);
			}

			std::vector<std::vector<int>>nodes_of_representative(interior_node_count);
			for(int x=0; x<interior_node_count; ++x)
				if(!is_in_sub_separator(x))
					nodes_of_representative[uf(x)].push_back(x);

			auto&separator = sub_separator;
			for(auto&x:separator)
				x = sub_node_to_node(x);

			for(int x=0; x<interior_node_count; ++x){
				if(!nodes_of_representative[x].empty()){
					Cell new_cell;

					auto&new_cell_interior_node_list = nodes_of_representative[x];
					for(auto&x:new_cell_interior_node_list)
						x = sub_node_to_node(x);

					new_cell.parent_cell = closed_cell_id;
					
					new_cell.separator_node_list = std::move(new_cell_interior_node_list);

					new_cell.boundary_node_list = current_cell.boundary_node_list;
					new_cell.boundary_node_list.insert(new_cell.boundary_node_list.end(), separator.begin(), separator.end());

					{
						for(auto x:new_cell.separator_node_list)
							in_child_cell.set(x, true);
						new_cell.boundary_node_list.erase(
							std::remove_if(
								new_cell.boundary_node_list.begin(),
								new_cell.boundary_node_list.end(),
								[&](int x)->bool{
									for(auto xy:inv_tail(x))
										if(in_child_cell(head(xy)))
											return false;
									return true;
								}
							),
							new_cell.boundary_node_list.end()
						);
						for(auto x:new_cell.separator_node_list)
							in_child_cell.set(x, false);
					}

					new_cell.separator_node_list.shrink_to_fit();
					new_cell.boundary_node_list.shrink_to_fit();
	
					if(new_cell.bag_size() > max_open_bag_size)
						max_open_bag_size = new_cell.bag_size();

					open_cells.push(std::move(new_cell));
				}
			}	

			current_cell.separator_node_list = std::move(separator);
			current_cell.separator_node_list.shrink_to_fit();

			for(int x:interior_node_list)
				node_to_sub_node[x] = -1;
		}

		if(current_cell.bag_size() > max_closed_bag_size)
			max_closed_bag_size = current_cell.bag_size();
		
		if(must_recompute_max_open_bag_size){
			max_open_bag_size = 0;
			for(auto&x:access_internal_vector(open_cells))
				if(x.bag_size() > max_open_bag_size)
					max_open_bag_size = x.bag_size();
		}

		closed_cells.push_back(std::move(current_cell));

		check_if_better();

		if(max_closed_bag_size >= smallest_known_treewidth){
			return;
		}

		if(max_closed_bag_size >= max_open_bag_size){
			return;
		}
	}
}

ArrayIDIDFunc preorder, inv_preorder;

std::string format_multilevel_partition_as_tree_decomposition(const std::vector<Cell>&cell_list){
	std::ostringstream out;
	print_tree_decompostion_of_multilevel_partition(out, tail, head, preorder, cell_list);
	return out.str();
}

void run_flow_cutter(int nb_nodes, const std::vector<int>& start_edges, const std::vector<int>& end_edges, 
	std::vector<std::vector<int>>& tree_bags,
	std::vector<int>& start_tree_edges,
	std::vector<int>& end_tree_edges,
	std::string method,
	int nb_runs, int random_seed){

	try{
		{
			auto g = load_graph(nb_nodes, start_edges, end_edges);
			tail = std::move(g.tail);
			head = std::move(g.head);
		}

		{
			preorder = compute_preorder(compute_successor_function(tail, head));
			for(int i=0; i<tail.image_count(); ++i)
				preorder[i] = i;
			inv_preorder = inverse_permutation(preorder);
			tail = chain(std::move(tail), inv_preorder);
			head = chain(std::move(head), inv_preorder);
		}

		{
			auto p = sort_arcs_first_by_tail_second_by_head(tail, head);
			tail = chain(p, std::move(tail));
			head = chain(p, std::move(head));
		}

		{
			try{
				std::minstd_rand rand_gen;
				rand_gen.seed(random_seed);
				ArrayIDIDFunc order;
				if (method == "min_shortcut"){
					order = chain(compute_greedy_min_shortcut_order(tail, head), inv_preorder);
				}
				else {
					order = chain(compute_greedy_min_degree_order(tail, head), inv_preorder);
				}
				store_tree_decompostion_of_order(tail, head, order, tree_bags, start_tree_edges, end_tree_edges);
			}catch(...){
			}
			
		}
	}catch(...){
	}
}

