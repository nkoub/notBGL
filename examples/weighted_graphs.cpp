/*
 *
 * This is an example of how to use the functions related to weighted graphs.
 *
 * Compilation requires the c++11 standard.
 *   Example: g++ -O3 -std=c++1y weighted_graphs.cpp
 * 
 * Author:  Antoine Allard
 * WWW:     antoineallard.info
 * Date:    May 2016
 * 
 */


// Standard template library 
#include <iostream>
// notBGL
#include "../src/notBGL.hpp"


int main(int argc, char** argv)
{
  
  // Instantiates an empty undirected graph.
  typedef boost::adjacency_list< boost::setS,
                                 boost::vecS,
                                 boost::undirectedS,
                                 notBGL::MinimalWeightedVertexProp,
                                 boost::property<boost::edge_weight_t, double> > graph_t;
  graph_t graph;
 
  // Populates the graph via the edgelist TestGraph2.edge.
  notBGL::load_weighted_edgelist("edgelists/TestGraph2.edge", graph);

  // Saves the weighted edgelist (see documentation of save_edges_properties() for further options).
  notBGL::save_weighted_edgelist("TestGraph2.edge", graph);

  // Extracts the strength of the vertices.
  auto Vertex2Strength = notBGL::strengths(graph);

  // Computes the disparity of vertices.
  auto Vertex2Disparity = notBGL::disparity(graph);

  // Writes the vertices' properties into a file (see documentation of save_vertices_properties() for further options)
  notBGL::save_vertices_properties("TestGraph2_vertices_prop.dat", graph, {Vertex2Strength, Vertex2Disparity});

  // Exits the program successfully.
  return 0;

  
}
