/*
 *
 * This is an example of how to use the functions related to clustering.
 *
 * Compilation requires the c++11 standard.
 *   Example: g++ -O3 -std=c++1y clustering.cpp
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
                                 notBGL::MinimalVertexProp> graph_t;
  graph_t graph;
 
  // Populates the graph via the edgelist TestGraph1.edge.
  notBGL::load_edgelist("edgelists/TestGraph1.edge", graph);

  // Surveys the triangles present in the graph. This function yields a std::vector whose elements
  // are std::tuple containing the boost::vertex_descriptor of the vertices participating to the
  // triangle.
  auto triangles = notBGL::survey_triangles(graph);

  // Computes the local clustering coefficients. Having identified the triangles in the graph, this
  // function computes the local clustering coefficient of each vertex and returns a std::map
  // mapping the boost::vertex_descriptor to the value of the coefficient.
  auto Vertex2LCC = notBGL::local_clustering_coefficients(triangles, graph);

  // Computes the multiplicity of edges. Having identified the triangles in the graph, this
  // function computes the multiplicity of each edge and returns a std::map
  // mapping the boost::edge_descriptor to the value of the multiplicity.
  auto Edge2Multiplicity = notBGL::multiplicity(triangles, graph);

  // Computes the global clustering coefficient, and prints it on screen.
  std::cout << "Global clustering coefficient: "
            << notBGL::global_clustering_coefficient(triangles, graph)
            << std::endl << std::endl;

  // Computes the average local clustering coefficient, and prints it on screen.
  std::cout << "Average local clustering coefficient: "
            << notBGL::average_local_clustering_coefficient(Vertex2LCC)
            << std::endl << std::endl;

  // Prints the local clustering coefficients on screen.
  std::cout << "Local clustering coefficients:" << std::endl;
  for(auto el : Vertex2LCC)
  {
    std::cout << graph[el.first].name << ": " << el.second << std::endl;
  }

  // Prints the multiplicities on screen.
  std::cout << std::endl << "Edge multiplicities:" << std::endl;
  for(auto el : Edge2Multiplicity)
  {
    std::cout << graph[boost::source(el.first, graph)].name << "  "
              << graph[boost::target(el.first, graph)].name << ": "
              << el.second << std::endl;
  }

  // Exits the program successfully.
  return 0;

  /* 
   * The output should read:
   * 
   * >> g++ -O3 -std=c++1y clustering.cpp -o clustering
   * >> ./clustering
   * Global clustering coefficient: 0.363636
   * 
   * Average local clustering coefficient: 0.539683
   * 
   * Local clustering coefficients:
   * Mark: 0.190476
   * Anna: 1
   * Tony: 0.666667
   * Nick: 1
   * Mary: 0
   * Lucy: 1
   * Fred: 0.666667
   * Vera: 0.333333
   * Nora: 0
   * 
   * Edge multiplicities:
   * Mark  Anna: 1
   * Mark  Tony: 2
   * Mark  Lucy: 1
   * Mark  Nick: 1
   * Mark  Mary: 0
   * Vera  Nora: 0
   * Mark  Fred: 2
   * Mark  Vera: 1
   * Anna  Tony: 1
   * Tony  Nick: 1
   * Lucy  Fred: 1
   * Fred  Vera: 1
   *
   */
}
