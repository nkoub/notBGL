/*
 *
 * This is an example of how to use the betweenness_centrality() function. The
 * code instantiates a simple and small undirected graph, computes the
 * betweenness centrality of its vertices and edges, and then prints the results
 * on the screen.
 *
 * Compilation requires the c++11 standard.
 *   Example: g++ -O3 -std=c++1y betweenness_centrality.cpp
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
  
  // Instantiates an undirected graph with 9 vertices.
  typedef boost::adjacency_list< boost::setS,
                                 boost::vecS,
                                 boost::undirectedS > graph_t;
  graph_t graph(9);
 
  // Creates the edges. Since we use boost::vecS to store the vertices, we can
  // directly use numerical IDs to identify each vertex.
  boost::add_edge(0, 1, graph);
  boost::add_edge(0, 2, graph);
  boost::add_edge(0, 3, graph);
  boost::add_edge(0, 4, graph);
  boost::add_edge(0, 5, graph);
  boost::add_edge(0, 6, graph);
  boost::add_edge(0, 7, graph);
  boost::add_edge(1, 2, graph);
  boost::add_edge(2, 3, graph);
  boost::add_edge(5, 6, graph);
  boost::add_edge(6, 7, graph);
  boost::add_edge(7, 8, graph);

  // Computes the betweenness centrality of the vertices and the edges. The
  // function returns a std::tuple of two std::map objects that map a
  // boost::vertex_descriptor and a boost::edge_descriptor to their
  // betweenness_centrlity, respectively.
  auto bc = notBGL::betweenness_centrality(graph);

  // Prints the betweenness centrality of vertices. We use std::get to access
  // the different map objects.
  std::cout << "Betweenness centrality of vertices" << std::endl;
  for(auto el : std::get<0>(bc))
  {
    std::cout << el.first << ": " << el.second << std::endl;
  }

  // Prints the betweenness centrality of edges.
  std::cout << std::endl << "Betweenness centrality of edges" << std::endl;
  for(auto el : std::get<1>(bc))
  {
    std::cout << boost::source(el.first, graph)
              << " <--> " 
              << boost::target(el.first, graph) 
              << ": " 
              << el.second 
              << std::endl;
  }

  // Exits the program successfully.
  return 0;

  /* 
   * The output should read:
   * 
   * >> g++ -O3 -std=c++1y betweenness_centrality.cpp -o betweenness_centrality
   * >> ./betweenness_centrality
   * Betweenness centrality of vertices
   * 0: 20.5
   * 1: 0
   * 2: 0.5
   * 3: 0
   * 4: 0
   * 5: 0
   * 6: 1
   * 7: 7
   * 8: 0
   *
   * Betweenness centrality of edges
   * 0 <--> 1: 6.5
   * 0 <--> 2: 6
   * 0 <--> 3: 6.5
   * 0 <--> 4: 8
   * 0 <--> 5: 6
   * 0 <--> 6: 5
   * 0 <--> 7: 11
   * 1 <--> 2: 1.5
   * 2 <--> 3: 1.5
   * 5 <--> 6: 2
   * 6 <--> 7: 3
   * 7 <--> 8: 8
   *
   */

}
