/*
 *
 * This is an example of how to use the functions related to finding connected components in undirected graphs.
 *
 * Compilation requires the c++11 standard.
 *   Example: g++ -O3 -std=c++1y connected_components.cpp
 * 
 * Author:  Antoine Allard
 * WWW:     antoineallard.info
 * Date:    July 2016
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
  notBGL::load_edgelist("edgelists/TestGraph4.edge", graph);

  // Identifies the connected components. The function returns a std::tuple<std::map<vertex_t>,
  // int>. The first element corresponds to a std::map<vertex_t, int> mapping each vertex to the
  // component it belongs to. The second elements is the number of components.
  auto comp = notBGL::connected_components(graph);

  // Computes the size of each component.
  unsigned int nb_components = std::get<1>(comp);
  std::vector<unsigned int> sizes(nb_components, 0);
  for(auto el : std::get<0>(comp))
  {
    sizes[el.second] += 1;
  }

  // Gather the vertices into the different components.
  unsigned int c;
  typedef typename graph_t::vertex_descriptor vertex_t;
  std::vector< std::vector<vertex_t> > components(nb_components);
  for(unsigned int c(0); c<nb_components; ++c)
  {
    components[c].resize(sizes[c]);
  }
  for(auto el : std::get<0>(comp))
  {
    c = el.second;
    sizes[c] -= 1;
    components[c][sizes[c]] = el.first;
  }


  // Prints the number of components on screen.
  std::cout << "Number of components: " << nb_components << std::endl << std::endl;

  // Prints the size of each component and its members.
  for(unsigned int c(0); c<nb_components; ++c)
  {
    std::cout << "Component " << c << " contains " << components[c].size()
              << " vertex/vertices." << std::endl;
    for(unsigned int v(0), vv(components[c].size()); v<vv; ++v)
    {
      std::cout << graph[components[c][v]].name << "  ";
    }
    std::cout << std::endl << std::endl;
  }


  // Exits the program successfully.
  return 0;

  /* 
   * The output should read:
   * 
   * >> g++ -O3 -std=c++1y connected_components.cpp -o connected_components
   * >> ./connected_components
   * Number of components: 3
   *
   * Component 0 contains 9 vertex/vertices.
   * Nora  Vera  Fred  Lucy  Mary  Nick  Tony  Anna  Mark  
   *
   * Component 1 contains 2 vertex/vertices.
   * Alba  Abby  
   *
   * Component 2 contains 5 vertex/vertices.
   * Isla  Fern  Erin  Carl  Clem
   *
   */
}
