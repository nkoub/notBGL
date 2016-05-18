/**
* @brief  Computes the absolute betweenness centrality of the vertices and the edges of a graph.
* @author Antoine Allard (<a href="http://antoineallard.info">antoineallard.info</a>)
* @date   May 2016
* @bug    No known bugs.
* @todo   Complete the documentation using doxygen.
*
*/
 
#include <iostream>
#include "../src/notBGL.hpp"
 
int main(int argc, char** argv)
{
  
  // Create the graph object.
  typedef boost::adjacency_list< boost::setS,
                                 boost::vecS,
                                 boost::undirectedS > graph_t;
  graph_t graph(9);
 
  // Populates the graph.
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

  // Computes the betweenness centrality of the vertices and the edges.
  typedef graph_t::vertex_descriptor vertex_t;
  typedef graph_t::edge_descriptor edge_t;
  std::map<vertex_t, double> Vertex2BC;
  std::map<edge_t, double> Edge2BC;
  std::tie(Vertex2BC, Edge2BC) = notBGL::betweenness_centrality(graph);

  // Prints the betweenness centrality of vertices.
  for(auto el : Vertex2BC)
  {
    std::cout << el.first << ": " << el.second << std::endl;
  }

  // Prints the betweenness centrality of edges.
  for(auto el : Edge2BC)
  {
    std::cout << boost::source(el.first, graph)
              << " <--> " 
              << boost::target(el.first, graph) 
              << ": " 
              << el.second 
              << std::endl;
  }

  return 0;
}
