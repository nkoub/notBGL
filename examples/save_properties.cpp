/*
 *
 * This is an example of how to use the function save_vertices_properties() function.
 *
 * Compilation requires the c++11 standard.
 *   Example: g++ -O3 -std=c++1y save_vertices_properties.cpp
 * 
 * Author:  Antoine Allard
 * WWW:     antoineallard.info
 * Date:    May 2016
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
  notBGL::load_edgelist("edgelists/TestGraph1.edge", graph);

  // Extracts the degree of vertices.
  auto Vertex2Degree = notBGL::degrees(graph);

  // Computes the betweenness centrality of the vertices and the edges. This function returns a
  // std::tuple of two std::map objects that map a boost::vertex_descriptor and a
  // boost::edge_descriptor to their betweenness_centrlity, respectively.
  auto BC = notBGL::betweenness_centrality(graph);

  // Surveys the triangles present in the graph. This function yields a std::vector whose elements
  // are std::tuple containing the boost::vertex_descriptor of the vertices participating to the
  // triangle.
  auto triangles = notBGL::survey_triangles(graph);
  // Computes the local clustering coefficients. Having identified the triangles in the graph, this
  // function computes the local clustering coefficient of each vertex and returns a std::map
  // mapping the boost::vertex_descriptor to the value of the coefficient.
  auto Vertex2LCC = notBGL::local_clustering_coefficients(triangles, graph);

  // Saves the properties of the vertices into a file. Note that the std::map related to the
  // properties to be written into the file are passed to the function as a std::initializer_list
  // object. 
  // By default, the name of the vertices are printed and that the width of the columns
  // are set to 10.
  notBGL::save_vertices_properties("TestGraph1_vertices_prop_name_width10_default.dat", graph, {Vertex2Degree, std::get<0>(BC), Vertex2LCC});
  // Identifies the vertices with a contiguous sequence of integers starting at 0.
  notBGL::save_vertices_properties("TestGraph1_vertices_prop_num_width10.dat", graph, {Vertex2Degree, std::get<0>(BC), Vertex2LCC}, notBGL::vertexID_num);
  // Only prints the properties without identifying the vertices, and changes the width of the
  // column to 15.
  notBGL::save_vertices_properties("TestGraph1_vertices_prop_none_width15.dat", graph, {Vertex2Degree, std::get<0>(BC), Vertex2LCC}, notBGL::vertexID_none, 15);

  // Same sequence of options, but for the edge properties.
  notBGL::save_edges_properties("TestGraph1_edges_prop_name_width10_default.dat", graph, {std::get<1>(BC)});
  notBGL::save_edges_properties("TestGraph1_edges_prop_num_width10.dat", graph, {std::get<1>(BC)}, notBGL::vertexID_num);
  notBGL::save_edges_properties("TestGraph1_edges_prop_none_width15.dat", graph, {std::get<1>(BC)}, notBGL::vertexID_none, 15);


  // Exits the program successfully.
  return 0;
}
