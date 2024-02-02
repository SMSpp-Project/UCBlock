#ifndef LIST_GRAPH_H
#define LIST_GRAPH_H

#include "array_id_func.h"

#include <tuple>
#include <vector>
#include <string>

struct ListGraph{
	ListGraph()=default;
	ListGraph(int node_count, int arc_count)
		:head(arc_count, node_count), tail(arc_count, node_count){}

	int node_count()const{ return head.image_count(); }
	int arc_count()const{ return head.preimage_count(); }

	ArrayIDIDFunc head, tail;
};

ListGraph load_graph(int node_count, const std::vector<int>& start_edges, const std::vector<int>& end_edges);

#endif
