/*
 *
 * This is an example of how to use the kcore_decomposition() and
 * mcore_decomposition() functions.
 *
 * Compilation requires the c++11 standard.
 *   Example: g++ -O3 -std=c++1y core_decomposition.cpp
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
  auto Vertex2kCore = notBGL::kcore_decomposition(graph);

  // Extract the m-coreness of vertices.
  auto triangles = notBGL::survey_triangles(graph);
  auto Edge2Multiplicity = notBGL::multiplicity(triangles, graph);
  auto Vertex2mCore = notBGL::mcore_decomposition(Edge2Multiplicity, graph);

  // Prints the corenesses of vertices.
  double width = 15;
  std::cout << std::setw(1) << "#";
  std::cout << std::setw(width-1) << "Vertex"     << " ";
  std::cout << std::setw(width)   << "k-coreness" << " ";
  std::cout << std::setw(width)   << "m-coreness" << " ";
  std::cout << std::endl;
  typename graph_t::vertex_iterator v_it, v_end;
  for(std::tie(v_it, v_end) = boost::vertices(graph); v_it!=v_end; ++v_it)
  {
    // Name.
    std::cout << std::setw(width) << graph[*v_it].name << " ";
    // k-core.
    std::cout << std::setw(width) << Vertex2kCore[*v_it] << " ";
    // m-core.
    std::cout << std::setw(width) << Vertex2mCore[*v_it] << " " << std::endl;
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
