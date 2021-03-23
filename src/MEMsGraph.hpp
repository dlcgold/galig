#ifndef _MEMSGRAPH_HPP_
#define _MEMSGRAPH_HPP_

#include <string>
#include <list>
#include <forward_list>
#include <utility>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <unordered_set>

#include "utils.hpp"
#include "SplicingGraph.hpp"
#include "bMEM.hpp"

#include <lemon/concepts/maps.h>
#include <lemon/list_graph.h>
#include <lemon/dijkstra.h>
#include <lemon/fib_heap.h>

typedef lemon::ListDigraph Graph;
typedef Graph::Node Node;
typedef Graph::NodeIt NodeIt;
typedef Graph::Arc Arc;
typedef Graph::ArcIt ArcIt;
typedef Graph::InArcIt InArc;
typedef Graph::OutArcIt OutArc;
typedef lemon::Path<Graph> Path;
typedef Graph::NodeMap<Mem> Node2MEM;
typedef Graph::ArcMap<int> Arc2Int;
typedef Graph::NodeMap<int> FibM;
typedef lemon::FibHeap<int, FibM> FibH;

class MemsGraph {
private:
  /**
   * L: L parameter 
   * eps: epsilon parameter
   * read: read
   * m: read length
   * K0: start/end proximity
   * K1: gap length
   * K2: errors in gap
   * exsN: exons number
   * graph: MEMsGraph
   * nodes_map: node -> mem
   * edges_map: edge -> weigth (int)
   * start: global starting node
   * end: global ending node
   * starting_nodes: local starting nodes list (exon)
   * ending_nodes: local ending nodes list (exon)
   **/
  int L; // lunghezza minima mem
  int LI;
  int eps;
  std::string read;
  bool verbose;
  bool intron_calc;
  int m;
  int K0;
  int K1;
  int K2;
  int exsN;
  Graph AnnGraph;
  Graph NovGraph;
  Node2MEM AnnNodesMap;
  Arc2Int AnnEdgesMap;
  Node2MEM NovNodesMap;
  Arc2Int NovEdgesMap;
  Node AnnStart;
  Node AnnEnd;
  Node NovStart;
  Node NovEnd;
  std::pair<bool, int> checkMEMs(const SplicingGraph&, const Mem&, const Mem&);
  std::pair<bool, int> validStart(const SplicingGraph&, const Mem&);
  std::pair<bool, int> validEnd(const SplicingGraph&, const Mem&);
  std::pair<bool, int> checkMEMsIntrons(const SplicingGraph&, const Mem&, const Mem&);
  std::pair<bool, int> validStartIntrons(const SplicingGraph&, const Mem&);
  std::pair<bool, int> validEndIntrons(const SplicingGraph&, const Mem&);

public: 

  MemsGraph(const std::string&,
	    const int&,
	    const int&,
	    const int&,
	    const bool&,
	    const bool&,
	    const bool&);
  void build(const SplicingGraph&,
	     const std::list<Mem>&,
	     BackwardMEM&);
  std::list<std::pair<int, std::list<Mem> > > visit(const SplicingGraph&);
  void save(const std::string&);
  
  std::vector<std::vector<Mem> > startIntrons(const Mem&, const SplicingGraph&,
					      std::list<Mem>&, const int&);
  std::vector<std::vector<Mem> > middleIntrons(const Mem&, const Mem&,
					       const SplicingGraph&,
					       std::list<Mem>&, const int&);
  std::vector<std::vector<Mem> > endIntrons(const Mem&, const SplicingGraph&,
					    std::list<Mem>&, const int&);
  std::vector<std::vector<Mem> > totIntrons(std::list<Mem>&);
};

#endif
