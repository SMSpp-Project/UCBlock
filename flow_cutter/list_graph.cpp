#include "list_graph.h"
#include "io_helper.h"
#include "multi_arc.h"
#include "id_multi_func.h"

#include <stdexcept>
#include <fstream>
#include <iostream>
#include <sstream>

ListGraph load_graph(int node_count, const std::vector<int>& start_edges, const std::vector<int>& end_edges){
	ListGraph graph;
	int arc_count = start_edges.size();
	graph = ListGraph(node_count, 2*arc_count);
	int next_arc = 0;

	std::vector<int>::const_iterator start_it = start_edges.begin();
	std::vector<int>::const_iterator end_it = end_edges.begin();
	while(start_it != start_edges.end() && end_it != end_edges.end()){
		int h = *start_it;
		int t = *end_it;
		if(h < 0 || h >= graph.node_count() || t < 0 || t >= graph.node_count())
			throw std::runtime_error("Invalid arc in line num "+std::to_string(next_arc));
		if(next_arc < graph.arc_count()){
			graph.head[next_arc] = h;
			graph.tail[next_arc] = t;
		}
		++next_arc;
		if(next_arc < graph.arc_count()){
			graph.head[next_arc] = t;
			graph.tail[next_arc] = h;
		}
		++next_arc;
		++start_it;
		++end_it;
	}

	if(next_arc != graph.arc_count())
		throw std::runtime_error("The arc count in the header ("+std::to_string(graph.arc_count())+") does not correspond with the actual number of arcs ("+std::to_string(next_arc)+").");

	return graph; // NVRO
}
/*
ListGraph uncached_load_pace_graph(const std::string&file_name){
	return load_uncached_text_file(file_name, load_pace_graph_impl);
}
*/

