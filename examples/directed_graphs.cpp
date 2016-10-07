/*
 *
 * This is an example of how to use the functions related to directed graphs.
 *
 * Compilation requires the c++11 standard.
 *   Example: g++ -O3 -std=c++1y directed_graphs.cpp
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
  
  // Instantiates an empty directed graph (bidirectional).
  typedef boost::adjacency_list< boost::setS,
                                 boost::vecS,
                                 boost::bidirectionalS,
                                 notBGL::MinimalVertexProp> graph_t;
  graph_t graph;
 
  // Populates the graph via the edgelist TestGraph3.edge.
  notBGL::load_edgelist("edgelists/TestGraph3.edge", graph);
  // notBGL::load_edgelist("edgelists/TestGraph5.edge", graph);

  // Extracts the out-degree of vertices.
  auto Vertex2OutDegree = notBGL::out_degrees(graph);

  // Extracts the in-degree of vertices.
  auto Vertex2InDegree = notBGL::in_degrees(graph);

  // Computes the number of reciprocal edge pairs for each vertex.
  auto Vertex2ReciprocalPairs = notBGL::reciprocical_edge_pairs(graph);

  // Computes the reciprocity. The function returns a std::tuple<double,double>. The first value
  // corresponds to the "naive definition" of reciprocity, that is the fraction of directed edges
  // that have a "twin" edge in the other direction (note that there is an even number of such
  // edges). The second coefficient is the coefficient defined in doi:10.1103/PhysRevLett.93.268701
  // that accounts for the reciprocity can occur randomly.
  auto Reciprocity = notBGL::reciprocity(Vertex2ReciprocalPairs, graph);

  // Extracts the triangles in the graph regardless of the nature of the connections.
  auto triangles = notBGL::survey_directed_triangles(graph);

  // Extracts the triangle spectrum.
  auto spectrum = notBGL::triangle_spectrum(triangles, graph);

  // Prints the reciprocity coefficients on screen..
  std::cout << "Naive reciprocity coefficient: "
            << std::get<0>(Reciprocity)
            << std::endl;
  std::cout << "Corrected reciprocity coefficient: "
            << std::get<1>(Reciprocity)
            << std::endl
            << std::endl;

  // Prints the degrees of nodes.
  double width = 15;
  std::cout << std::setw(1) << "#";
  std::cout << std::setw(width-1) << "Vertex"        << " ";
  std::cout << std::setw(width)   << "In-degree"     << " ";
  std::cout << std::setw(width)   << "Out-degree"    << " ";
  std::cout << std::setw(width)   << "Recip. pairs"  << " ";
  std::cout << std::setw(width)   << "Undir. degree" << " ";
  std::cout << std::endl;
  typename graph_t::vertex_iterator v_it, v_end;
  for(std::tie(v_it, v_end) = boost::vertices(graph); v_it!=v_end; ++v_it)
  {
    // Name.
    std::cout << std::setw(width) << graph[*v_it].name << " ";
    // In-degree.
    std::cout << std::setw(width) << Vertex2InDegree[*v_it] << " ";
    // Out-degree.
    std::cout << std::setw(width) << Vertex2OutDegree[*v_it] << " ";
    // Number of reciprocal pairs.
    std::cout << std::setw(width) << Vertex2ReciprocalPairs[*v_it] << " ";
    // Degree in the undirected projection of the directed graph.
    std::cout << std::setw(width) << Vertex2InDegree[*v_it] +
                                     Vertex2OutDegree[*v_it] -
                                     Vertex2ReciprocalPairs[*v_it] << " ";
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // Prints the triangles in the graph.
  bool is_edge;
  typename graph_t::edge_descriptor e;
  std::cout << "Triangles:" << std::endl;
  for(auto t : triangles)
  {
    std::cout << graph[t[0]].name;
    std::tie(e, is_edge) = boost::edge(t[1], t[0], graph);
    is_edge ? std::cout << " <" : std::cout << " -";
    std::cout << "--";
    std::tie(e, is_edge) = boost::edge(t[0], t[1], graph);
    is_edge ? std::cout << "> " : std::cout << "- ";
    std::cout << graph[t[1]].name;
    std::tie(e, is_edge) = boost::edge(t[2], t[1], graph);
    is_edge ? std::cout << " <" : std::cout << " -";
    std::cout << "--";
    std::tie(e, is_edge) = boost::edge(t[1], t[2], graph);
    is_edge ? std::cout << "> " : std::cout << "- ";
    std::cout << graph[t[2]].name;
    std::tie(e, is_edge) = boost::edge(t[0], t[2], graph);
    is_edge ? std::cout << " <" : std::cout << " -";
    std::cout << "--";
    std::tie(e, is_edge) = boost::edge(t[2], t[0], graph);
    is_edge ? std::cout << "> " : std::cout << "- ";
    std::cout << graph[t[0]].name << std::endl;
  }
  std::cout << std::endl;

  // Prints the histogram of triangle types.
  std::cout << "Number of each triangle types:" << std::endl;
  std::cout << "A ---> B ---> C ---> A: "<< spectrum[0] << std::endl;
  std::cout << "A ---> B ---> C <--- A: "<< spectrum[1] << std::endl;
  std::cout << "A ---> B <--- C <--> A: "<< spectrum[2] << std::endl;
  std::cout << "A ---> B ---> C <--> A: "<< spectrum[3] << std::endl;
  std::cout << "A <--- B ---> C <--> A: "<< spectrum[4] << std::endl;
  std::cout << "A ---> B <--> C <--> A: "<< spectrum[5] << std::endl;
  std::cout << "A <--> B <--> C <--> A: "<< spectrum[6] << std::endl;



  // Exits the program successfully.
  return 0;

  /* 
   * The output should read:
   * 
   * >> g++ -O3 -std=c++1y directed_graph.cpp -o directed_graph
   * >> ./directed_graph
   * Naive reciprocity coefficient: 0.666667
   * Corrected reciprocity coefficient: 0.555556
   *
   * #        Vertex       In-degree      Out-degree    Recip. pairs   Undir. degree 
   *            Mark               3               6               2               7 
   *            Anna               2               2               2               2 
   *            Tony               2               3               2               3 
   *            Nick               2               0               0               2 
   *            Lucy               2               1               1               2 
   *            Fred               3               2               2               3 
   *            Vera               3               2               2               3 
   *            Mary               0               1               0               1 
   *            Nora               1               1               1               1
   *
   * Triangles:
   * Mark <-> Anna <-> Tony <-> Mark
   * Mark <-> Tony --> Nick <-- Mark
   * Mark --> Lucy <-> Fred <-- Mark
   * Mark --> Fred <-> Vera <-- Mark
   *
   * Number of each triangle types:
   * A ---> B ---> C ---> A: 0
   * A ---> B ---> C <--- A: 0
   * A ---> B <--- C <--> A: 1
   * A ---> B ---> C <--> A: 0
   * A <--- B ---> C <--> A: 2
   * A ---> B <--> C <--> A: 0
   * A <--> B <--> C <--> A: 1
   *
   */
}
