/*
 *
 * This is an example of how to use the kcore_decomposition() function.
 *
 * Compilation requires the c++11 standard.
 *   Example: g++ -O3 -std=c++1y kcore_decomposition.cpp
 * 
 * Author:  Antoine Allard
 * WWW:     antoineallard.info
 * Date:    Oct. 2016
 * 
 */


// notBGL
#include "../src/notBGL.hpp"


int main(int argc, char** argv)
{
  
  // Instantiates an empty undirected graph.
  typedef boost::adjacency_list< boost::setS,
                                 boost::vecS,
                                 boost::undirectedS,
                                 notBGL::MinimalVertexProp> graph_t;
  graph_t graph;
 
  // Populates the graph via the edgelist TestGraph1.edge.
  notBGL::load_edgelist("edgelists/TestGraph6.edge", graph);

  // Extract the coreness of vertices.
  auto Vertex2Coreness = notBGL::kcore_decomposition(graph);

  // Prints the corenesses on screen.
  std::cout << "Vertex coreness:" << std::endl;
  for(auto el : Vertex2Coreness)
  {
    std::cout << graph[el.first].name << ": "
              << el.second << std::endl;
  }

  // Exits the program successfully.
  return 0;

  /* 
   * The output should read:
   * 
   * >> g++ -O3 -std=c++1y kcore_decomposition.cpp -o kcore_decomposition
   * >> ./kcore_decomposition
   * Vertex coreness:
   * Mark: 4
   * Anna: 4
   * Tony: 4
   * Nick: 4
   * Mary: 3
   * Lucy: 4
   * Vera: 2
   * Fred: 2
   * Nora: 1
   *
   */
}
