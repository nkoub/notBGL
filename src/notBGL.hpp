/**
 * @file       notBGL.hpp
 *
 * @brief      Source code of the notBGL library.
 *
 *             This file contains the complete source code of the notBGL
 *             library. While having one single file is not the most clear and
 *             organized choice for source codes, this format has been chosen to
 *             faciliate protability of the code.
 *
 * @author     Antoine Allard (<a
 *             href="http://antoineallard.info">antoineallard.info</a>)
 *
 * @date       March 2016
 */

#ifndef NOTBGL_HPP_INCLUDED
#define NOTBGL_HPP_INCLUDED

// Standard Template Library
#include <algorithm>         // std::transform
#include <cmath>             // std::sqrt, std::pow, std::cosh, std::sinh, std::acosh, std::fabs
#include <cstdlib>           // std::terminate
#include <fstream>           // std::ifstream, std::ofstream
#include <functional>        // std::minus
#include <initializer_list>  // std::initializer_list
#include <iomanip>           // std::setw
#include <iostream>          // std::cerr
#include <limits>            // std::numeric_limits
#include <list>              // std::list
#include <map>               // std::map
#include <numeric>           // std::inner_product
#include <queue>             // std::queue
#include <set>               // std::set
#include <sstream>           // std::stringstream
#include <string>            // std::string, std::stod, std::getline
#include <tuple>             // std::tie, std::make_tuple
#include <vector>            // std::vector
// Boost Graph Library
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/multi_array.hpp>

// All the methods are under the namespace "notBGL"
namespace notBGL
{

  #ifndef DOXYGEN_EXCLUDED
    // Minimal vertex property structure used in the generic undirected unweighted graph type (not documented).
    struct MinimalVertexProp  { std::string name; unsigned int num; };
    struct MinimalWeightedVertexProp  { std::string name; unsigned int num; double strength = 0; };

    // Default width of columns in output files.
    const unsigned int DEFAULT_COLUMN_WIDTH = 10;
  #endif // DOXYGEN_EXCLUDED

  /// Generic undirected unweighted graph type.
  typedef boost::adjacency_list<boost::listS,
                                boost::vecS,
                                boost::undirectedS,
                                notBGL::MinimalVertexProp,
                                boost::no_property > UndirectedGraph_t;

  /// Generic directed unweighted graph type.
  typedef boost::adjacency_list<boost::listS,
                                boost::vecS,
                                boost::bidirectionalS,
                                notBGL::MinimalVertexProp,
                                boost::no_property > DirectedGraph_t;

  /// Generic undirected weighted graph type.
  typedef boost::adjacency_list<boost::listS,
                                boost::vecS,
                                boost::undirectedS,
                                notBGL::MinimalWeightedVertexProp,
                                boost::property<boost::edge_weight_t, double> > WeightedGraph_t;

#ifndef DOXYGEN_EXCLUDED
  enum vertex_identifier_t {vertexID_none, vertexID_num, vertexID_name};

  
#endif // DOXYGEN_EXCLUDED

  /** ---------------------------------------------------------------------------------------------
   * @defgroup   IO Input / Output
   * @brief      This group contains the functions related to inputs and outputs
   *             such as loading graphs from edgelists from files as well as
   *             writing edgelists and proporties of vertices or edges into
   *             files.
   */

  /**
   * @brief      Populates a graph from an edge list.
   *
   * @param      filename  Path to the file containing the edgelist. The
   *                       edgelist file must have the format... anything at the
   *                       right will be ignored. In other words, a weighted
   *                       edgelist of the format... can be used without
   *                       problem. Lines beginning with "#" are ignored.
   * @param      g         Graph object to populate using an edgelist.
   *
   * @tparam     graph_t   { description }
   *
   * @return     Returns a std::map mapping the name of vertices to their
   *             vertex_descriptor. This object is provided to
   *
   * @ingroup    IO
   *
   * @see        save_edgelist()
   */
  template<typename graph_t>
  auto load_edgelist(std::string filename, graph_t &g);

  /**
   * @brief      Populates a graph from an edge list.
   *
   * @param      filename     Path to the file to write the edgelist into. The
   *                          edgelist file will have the format
   * @param      g            Graph object.
   * @param      write_names  If set to false, the numerical IDs of the vertices
   *                          are written instead of their names.
   *
   * @tparam     graph_t      { description }
   *
   * @ingroup    IO
   *
   * @see        load_edgelist()
   */
  template<typename graph_t>
  void save_edgelist(std::string filename, graph_t &g, bool write_names = true);

  /**
   * @brief      Saves properties of vertices into a file.
   *
   *             This is the main function to save properties of vertices into a
   *             file. Several other vertex-related output functions are
   *             wrappers of the save_vertices_properties(). See the
   *             `save_vertices_properties.cpp` example below for further
   *             details.
   *
   * @param      filename  Name of the file in which the properties are to
   *                       written.
   * @param      graph     The graph object.
   * @param      props     An std::initializer_list containing std::map objects
   *                       mapping boost::vertex_descriptor objects to the
   *                       values of the property of the vertices.
   * @param      vertexID  [optional] Identifier used for the vertices. Options
   *                       are
   *                       + `notBGL::vertexID_name` (default): prints the name
   *                         of the vertices;
   *                       + `notBGL::vertexID_num`: identifies the vertices
   *                         with a contiguous sequence of integers starting at
   *                         0;
   *                       + `notBGL::vertexID_none`: properties are printed
   *                         without identifying the vertices.
   * @param      width     [optional] Width of the columns in the file. Default
   *                       is 10.
   * @param      header    [optional] An std::initializer_list containing
   *                       std::string corresponding to the header of the column
   *                       of each property. A header "Vertex" is also added to
   *                       the column identifying the vertices, if applicable.
   *                       Also, a symbol "#" is added at the beginning of the
   *                       header to facilitate the use of the output file with
   *                       softwares like Python/Numpy. Default is no header.
   *
   * @tparam     graph_t   boost::adjacency_list
   * @tparam     vertex_t  boost::vertex_descriptor
   *
   * @note       For convenience, there exists wrappers of this function that
   *             allows the interchange of the parameters. See the
   *             `save_properties.cpp` example below and the notBGL.hpp source
   *             file for details.
   *
   * @see        save_edges_properties()
   *
   * @ingroup    IO
   */
  template<typename graph_t, typename vertex_t>
  void save_vertices_properties(std::string filename, graph_t &graph, std::initializer_list<std::map<vertex_t, double> > props, vertex_identifier_t vertexID = vertexID_name, unsigned int width = notBGL::DEFAULT_COLUMN_WIDTH, std::initializer_list<std::string> header = std::initializer_list<std::string>());

  template<typename graph_t, typename vertex_t>
  void save_vertices_properties(std::string filename, graph_t &graph, std::initializer_list<std::map<vertex_t, double> > props, std::initializer_list<std::string> header, vertex_identifier_t vertexID = vertexID_name, unsigned int width = notBGL::DEFAULT_COLUMN_WIDTH);

  /**
   * @brief      Saves properties of edges into a file.
   *
   *             This is the main function to save properties of edges into a
   *             file. Several other edge-related output functions are wrappers
   *             of the save_edges_properties(). See the `save_properties.cpp`
   *             example below for further details.
   *
   * @param      filename  Name of the file in which the properties are to
   *                       written.
   * @param      graph     The graph object.
   * @param      props     An std::initializer_list containing std::map objects
   *                       mapping boost::edge_descriptor objects to the values
   *                       of the property of the edges.
   * @param      vertexID  [optional] Identifier used for the edges. Options are
   *                       + `notBGL::vertexID_name` (default): prints the name
   *                         of the vertices;
   *                       + `notBGL::vertexID_num`: identifies the vertices
   *                         with a contiguous sequence of integers starting at
   *                         0;
   *                       + `notBGL::vertexID_none`: properties are printed
   *                         without identifying the edges.
   * @param      width     [optional] Width of the columns in the file. Default
   *                       is 10.
   * @param      header    [optional] An std::initializer_list containing
   *                       std::string corresponding to the header of the column
   *                       of each property. A header "Vertex" is also added to
   *                       the columns identifying the vertices, if applicable.
   *                       Also, a symbol "#" is added at the beginning of the
   *                       header to facilitate the use of the output file with
   *                       softwares like Python/Numpy. Default is no header.
   *
   * @tparam     graph_t   boost::adjacency_list
   * @tparam     edge_t    boost::edge_descriptor
   *
   * @note       For convenience, there exists wrappers of this function that
   *             allows the interchange of the parameters. See the
   *             `save_properties.cpp` example below and the notBGL.hpp source
   *             file for details.
   *
   * @see        save_vertices_properties()
   *
   * @ingroup    IO
   */
  template<typename graph_t, typename edge_t>
  void save_edges_properties(std::string filename, graph_t &graph, std::initializer_list<std::map<edge_t, double> > props, vertex_identifier_t vertexID = vertexID_name, unsigned int width = notBGL::DEFAULT_COLUMN_WIDTH, std::initializer_list<std::string> header = std::initializer_list<std::string>());


  template<typename graph_t, typename edge_t>
  void save_edges_properties(std::string filename, graph_t &graph, std::initializer_list<std::map<edge_t, double> > props, std::initializer_list<std::string> header, vertex_identifier_t vertexID = vertexID_name, unsigned int width = notBGL::DEFAULT_COLUMN_WIDTH);





  /** ---------------------------------------------------------------------------------------------
   * @defgroup   topo Topology
   * @brief      This group contains the functions related to the
   *             characterization of the topology of graphs (e.g., clustering,
   *             betweenness centrality).
   */

  /**
   * @brief      Extracts the degree sequence.
   *
   *             Extracts the degree sequence.
   *
   * @param      graph    The graph object
   *
   * @tparam     graph_t  boost::adjacency_list
   *
   * @return     Returns a std::map object mapping the boost:vertex_descriptor
   *             to the degree (double).
   *
   * @ingroup    topo
   */
  template<typename graph_t>
  auto degrees(graph_t &graph);

  /**
   * @brief      Identifies the triangles in the graph.
   *
   * @param      g        Graph object.
   *
   * @tparam     graph_t  { description }
   *
   * @return     A std::vector of tuple of vertex_descriptor
   *
   * @see        local_clustering_coefficient
   *
   * @ingroup    topo
   */
  template<typename graph_t>
  auto survey_triangles(graph_t &g);

  /**
   * @brief      Computes the local clustering coefficients.
   *
   *             Ceci est un `test`
   *
   * @param      triangles  Vector containing the triangles in the graph (from
   *                        survey_triangles()).
   * @param      graph      The graph
   *
   * @tparam     vector_t   { description }
   * @tparam     graph_t    { description }
   *
   * @return     The local clustering coefficient.
   *
   * @see        survey_triangles
   *
   * @ingroup    topo
   */
  template<typename vector_t, typename graph_t>
  auto local_clustering_coefficients(vector_t &triangles, graph_t &graph);

  /**
   * @brief      { function_description }
   *
   * @param[in]  local_clustering_coefficients  The local clustering
   *                                            coefficients
   *
   * @tparam     map_t                          { description }
   *
   * @return     { description_of_the_return_value }
   *
   * @ingroup    topo
   */
  template<typename map_t>
  double average_local_clustering_coefficient(map_t &local_clustering_coefficients);

  /**
   * @brief      Computes the global clustering coefficient.
   *
   *             Computes the global clustering coefficient according to
   *             @f[ C = \frac{3N_\triangle}{N_\wedge}
   *             @f] where
   *             @f$N_\triangle\f$ and
   *             @f$N_\wedge\f$ are the number of triangles and the number of
   *             triplets present in the graph, respectively. More details can
   *             be found
   *             [here](https://en.wikipedia.org/wiki/Clustering_coefficient#Global_clustering_coefficient)
   *
   * @param      triangles  Vector containing the triangles in the graph
   *                        [obtained from survey_triangles()].
   * @param      graph      The graph object.
   *
   * @tparam     vector_t   std::vector< std::tuple< boost::vertex_descriptor> > >
   * @tparam     graph_t    boost::adjacency_list
   *
   * @return     Value of the global clustering coefficient.
   *
   * @see        survey_triangles(), average_local_clustering_coefficient()
   * 
   * @ingroup    topo
   */
  template<typename vector_t, typename graph_t>
  double global_clustering_coefficient(vector_t &triangles, graph_t &graph);

  /**
   * @brief      Computes the multiplicity of edges.
   *
   *             The multiplicity of an edge is the number of triangles to which it participates.
   *
   * @param      triangles  Vector containing the triangles in the graph
   *                        [obtained from survey_triangles()].
   * @param      graph      The graph object.
   *
   * @tparam     vector_t   std::vector< std::tuple< boost::vertex_descriptor> >
   *                        >
   * @tparam     graph_t    boost::adjacency_list
   *
   * @return     std::map object mapping the boost::edge_descriptor to the
   *             multiplicities.
   *
   * @ingroup    topo
   * 
   * @see        survey_triangles()
   */
  template<typename vector_t, typename graph_t>
  auto multiplicity(vector_t &triangles, graph_t &graph);

  /**
   * @brief      Computes the betweenness centrality of the vertices and the
   *             edges of a graph.
   *
   * @param      graph    The graph object.
   *
   * @tparam     graph_t  boost::adjacency_list
   *
   * @return     { description_of_the_return_value }
   *
   * @ingroup    topo
   */
  template<typename graph_t>
  auto betweenness_centrality(graph_t &graph);

  /**
   * @brief      Identifies the components in which the vertices are.
   *
   * @param      graph    The graph
   *
   * @tparam     graph_t  { description }
   *
   * @return     { description_of_the_return_value }
   *
   * @ingroup    topo
   */
  template<typename graph_t>
  auto connected_components(graph_t &graph);


  /**
   * @brief      { function_description }
   *
   * @param      graph    The graph
   *
   * @tparam     graph_t  { description }
   *
   * @return     { description_of_the_return_value }
   * 
   * @ingroup    topo
   * 
   */
  template<typename graph_t>
  auto kcore_decomposition(graph_t &graph);



  /** ---------------------------------------------------------------------------------------------
   * @defgroup   directed Directed graphs
   * @brief      This group contains the functions specific to
   *             directed graphs (e.g., reciprocity, in/out-degree).
   */

   /**
   * @brief      Returns the in-degree sequence.
   *
   * @param      graph    The graph object
   *
   * @tparam     graph_t  boost::adjacency_list
   *
   * @return     { description_of_the_return_value }
   * 
   * @ingroup    directed
   */
  template<typename graph_t>
  auto in_degrees(graph_t &graph);

  /**
   * @brief      Returns the out-degree sequence.
   *
   * @param      graph    The graph object
   *
   * @tparam     graph_t  boost::adjacency_list
   *
   * @return     { description_of_the_return_value }
   * 
   * @ingroup    directed
   */
  template<typename graph_t>
  auto out_degrees(graph_t &graph);


  /**
   * @brief      { function_description }
   *
   * @param      graph    The graph object.
   *
   * @tparam     graph_t  { description }
   *
   * @return     { description_of_the_return_value }
   * 
   * @ingroup    directed
   */
  template<typename graph_t>
  auto reciprocical_edge_pairs(graph_t &graph);

  /**
   * @brief      { function_description }
   *
   * @param      Vector2ReciprocicalEdges  The vector 2 reciprocical edges
   * @param      graph                     The graph
   *
   * @tparam     map_t                     { description }
   * @tparam     graph_t                   { description }
   *
   * @return     { description_of_the_return_value }
   * 
   * @ingroup    directed
   */
  template<typename map_t, typename graph_t>
  auto reciprocity(map_t &Vector2ReciprocicalEdges, graph_t &graph);


  /**
   * @brief      { function_description }
   *
   * @param      g        { parameter_description }
   *
   * @tparam     graph_t  { description }
   *
   * @return     { description_of_the_return_value }
   */
  template<typename graph_t>
  auto survey_directed_triangles(graph_t &g);

  /**
   * @brief      { function_description }
   *
   * @param      triangles  The triangles
   * @param      graph      The graph
   *
   * @tparam     vector_t   { description }
   * @tparam     graph_t    { description }
   *
   * @return     { description_of_the_return_value }
   */
  template<typename vector_t, typename graph_t>
  std::vector<double> triangle_spectrum(vector_t &triangles, graph_t &graph);



  /** ---------------------------------------------------------------------------------------------
   * @defgroup   weights Weighted organization
   * @brief      This group contains the functions related to graphs embedded in
   *             a geometric space (e.g., loading coordinates, calculating
   *             distances, greedy routing).
   */

  /**
   * @brief      Extracts the strength sequence.
   *
   * @param      graph    The graph object.
   *
   * @tparam     graph_t  boost::adjacency_list
   *
   * @return     std::map object mapping the boost::vertex_descriptor to the
   *             strengths.
   *
   * @see        degrees()
   * 
   * @ingroup    weights
   */
  template<typename graph_t>
  auto strengths(graph_t &graph);

  /**
   * @brief      Extracts the weight sequence.
   *
   * @param      graph    The graph object.
   *
   * @tparam     graph_t  boost::adjacency_list
   *
   * @return     std::map object mapping the boost::edge_descriptor to the
   *             multiplicities.
   *
   * @ingroup    weights
   */
  template<typename graph_t>
  auto weights(graph_t &graph);

  /**
   * @brief      Computes the disparity of vertices.
   *
   * @param      graph    The graph object.
   *
   * @tparam     graph_t  boost::adjacency_list
   *
   * @return     std::map object mapping the boost::vertex_descriptor to the
   *             disparities.
   *
   * @ingroup    weights
   */
  template<typename graph_t>
  auto disparity(graph_t &graph);

  /**
   * @brief      Populates a graph from a weighted edge list.
   *
   * @param      filename  Path to the file containing the edgelist. The
   *                       edgelist file must have the format... anything at the
   *                       right will be ignored. In other words, a weighted
   *                       edgelist of the format... can be used without
   *                       problem.
   * @param      g         Graph object to populate using an edgelist.
   *
   * @tparam     graph_t   { description }
   *
   * @return     Returns a std::map mapping the name of vertices to their
   *             vertex_descriptor. This object is provided to
   *
   * @ingroup    weights
   *
   * @see        save_edgelist()
   */
  template<typename graph_t>
  auto load_weighted_edgelist(std::string filename, graph_t &g);

  /**
   * @brief      Saves the weighted edgelist into a file.
   *
   *             This function saves the weighted edgelist into a file. It is a
   *             wrapper of the more general function save_edges_properties().
   *             See the `save_properties.cpp` example below for further
   *             details.
   *
   * @param      filename  Name of the file in which the weighted edgelist is to
   *                       be written.
   * @param      graph     The graph object.
   * @param      vertexID  [optional] Identifier used for the edges. Options are
   *                       + `notBGL::vertexID_name` (default): prints the name
   *                         of the vertices;
   *                       + `notBGL::vertexID_num`: identifies the vertices
   *                         with a contiguous sequence of integers starting at
   *                         0;
   *                       + `notBGL::vertexID_none`: properties are printed
   *                         without identifying the edges.
   * @param      width     [optional] Width of the columns in the file. Default
   *                       is 10.
   *
   * @tparam     graph_t   boost::adjacency_list
   *
   * @see        save_edges_properties()
   *
   * @ingroup    weights
   */
  template<typename graph_t>
  void save_weighted_edgelist(std::string filename, graph_t &graph, vertex_identifier_t vertexID = vertexID_name, unsigned int width = 10);





  /** ---------------------------------------------------------------------------------------------
   * @defgroup   geo Geometry
   * @brief      This group contains the functions related to graphs embedded in
   *             a geometric space (e.g., loading coordinates, calculating
   *             distances, greedy routing).
   */

  /**
   * @brief      Loads the coordinates of vertices.
   *
   * @param[in]  filename     The file containing the coordinates in the format
   *                          [VertexName, x_1, x_2, x_3, ...]
   * @param      graph        The graph object.
   * @param      Name2Vertex  The name2 vertex
   *
   * @tparam     graph_t      { description }
   * @tparam     map_t        { description }
   *
   * @return     { description_of_the_return_value }
   *
   * @ingroup    geo
   */
  template<typename graph_t, typename map_t>
  auto load_coordinates(std::string filename, graph_t &graph, map_t &Name2Vertex);

  /**
   * @brief      { function_description }
   *
   * @param      x1    { parameter_description }
   * @param[in]  x2    { parameter_description }
   *
   * @return     { description_of_the_return_value }
   */
  double euclidean_distance(std::vector<double> &x1, std::vector<double> x2);

  /**
   * @brief      { function_description }
   *
   * @param      x1    { parameter_description }
   * @param[in]  x2    { parameter_description }
   * @param[in]  zeta  The zeta
   *
   * @return     { description_of_the_return_value }
   */
  double hyperbolic_distance(std::vector<double> &x1, std::vector<double> x2, double zeta);


  



  #ifndef DOXYGEN_EXCLUDED
    namespace utilities
    {
      template<typename map_t>
      double sum_of_map(map_t &map);

      template<typename map_t>
      double average_of_map(map_t &map);
    }
  #endif // DOXYGEN_EXCLUDED

}















// ================================================================================================
// ================================================================================================
// ================================================================================================
// *** Definition of functions.
// ================================================================================================










#ifndef DOXYGEN_EXCLUDED
// ================================================================================================
// ================================================================================================
template<typename map_t>
double notBGL::utilities::sum_of_map(map_t &map)
{
  double sum = 0;
  for(auto el : map)
  {
    sum += el.second;
  }
  return sum;
}





// ================================================================================================
// ================================================================================================
template<typename map_t>
double notBGL::utilities::average_of_map(map_t &map)
{
  // double sum = 0;
  // for(auto el : map)
  // {
  //   sum += el.second;
  // }
  // return sum / map.size();
  return notBGL::utilities::sum_of_map(map) / map.size();
}
#endif // DOXYGEN_EXCLUDED









// ================================================================================================
// ================================================================================================
// *** Module: Input/Output
// ================================================================================================


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::load_edgelist(std::string filename, graph_t &g)
{
  // Vertex descriptors.
  typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
  vertex_t v1, v2;
  // Stream objects.
  std::ifstream edgelist_file;
  std::stringstream one_line;
  // String objects.
  std::string full_line, name1_str, name2_str;
  // Integer objects.
  unsigned int node_cnt(0);
  // Map objects to assure the uniqueness of vertices.
  std::map< std::string, vertex_t > Name2Vertex;
  // Iterator objects.
  typename std::map< std::string, vertex_t >::iterator name_it;
  // Opens the stream and terminates if the operation did not succeed.
  edgelist_file.open(filename.c_str(), std::ios_base::in);
  if( !edgelist_file.is_open() )
  {
    std::cerr << "Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }
  else
  {
    // Reads the edgelist and registers the nodes and the edges.
    while( !edgelist_file.eof() )
    {
      // Reads a line of the edgelist.
      std::getline(edgelist_file,full_line); edgelist_file >> std::ws;
      one_line.str(full_line); one_line >> std::ws;
      one_line >> name1_str >> std::ws;
      // Skips a line of comment.
      if(name1_str == "#")
      {
        one_line.clear();
        continue;
      }
      one_line >> name2_str >> std::ws;
      one_line.clear();
      // Is name1 new?
      name_it = Name2Vertex.find(name1_str);
      if( name_it == Name2Vertex.end() )
      {
        v1 = boost::add_vertex(g);
        Name2Vertex[name1_str] = v1;
        g[v1].name = name1_str;
        g[v1].num = node_cnt;
        ++node_cnt;
      }
      else
      {
        v1 = name_it->second;
      }
      // Is name2 new?
      name_it = Name2Vertex.find(name2_str);
      if( name_it == Name2Vertex.end() )
      {
        v2 = boost::add_vertex(g);
        Name2Vertex[name2_str] = v2;
        g[v2].name = name2_str;
        g[v2].num = node_cnt;
        ++node_cnt;
      }
      else
      {
        v2 = name_it->second;
      }
      // Creates the edge if it does not exist (forbids multiple edges and self-loops).
      if(v1 != v2) // not a self-loop.
      {
        if(boost::edge(v1,v2,g).second == false) // not a multiple edge
        {
          boost::add_edge(v1, v2, g);
        }
      }
    }
  }
  // Closes the stream.
  edgelist_file.close();
  // Returns the maps "Name2Vertex" to help add other properties to vertices.
  return Name2Vertex;
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
void notBGL::save_edgelist(std::string filename, graph_t &g, bool write_names = true)
{
  // Stream objects.
  std::ofstream edgelist_file;
  // Opens the stream and terminates if the operation did not succeed.
  edgelist_file.open(filename.c_str(), std::ios_base::out);
  if( !edgelist_file.is_open() )
  {
    std::cerr << "Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }
  else
  {
    // Iterators over the edges of the graph.
    typename boost::graph_traits<graph_t>::edge_iterator it, end;
    // Writes each edge under the format source target, one edge per line.
    if(write_names)  // vertices are identified using their names.
    {
      for (std::tie(it, end) = edges(g); it!=end; ++it)
      {
        edgelist_file << g[source(*it, g)].name << " "
                      << g[target(*it, g)].name << std::endl;
      }
    }
    else // vertices are identified using the numerical id.
    {
      for (std::tie(it, end) = edges(g); it!=end; ++it)
      {
        edgelist_file << g[source(*it, g)].num << " "
                      << g[target(*it, g)].num << std::endl;
      }
    }
  }
  // Closes the stream.
  edgelist_file.close();
}


// ================================================================================================
// ================================================================================================
template<typename graph_t, typename vertex_t>
void notBGL::save_vertices_properties(std::string filename, graph_t &graph, std::initializer_list<std::map<vertex_t, double> > props, vertex_identifier_t vertexID, unsigned int width, std::initializer_list<std::string> header)
{
  // Stream objects.
  std::ofstream output_file;
  // Opens the stream and terminates if the operation did not succeed.
  output_file.open(filename.c_str(), std::ios_base::out);
  if( !output_file.is_open() )
  {
    std::cerr << "Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }
  else
  {
    // Iterators over the vertices of the graph.
    typename boost::graph_traits<graph_t>::vertex_iterator v_it, v_end;
    // Iterators over the std::string in "header".
    typename std::initializer_list<std::string>::iterator h_it(header.begin()), h_end(header.end());
    // Iterators over the std::map<vertex_t, double> in "props".
    typename std::initializer_list<std::map<vertex_t, double> >::iterator p_it, p_end(props.end());

    // Prints the vertex properties.
    if(vertexID == notBGL::vertexID_name)
    {
      if(header.size() > 0)
      {
        output_file << "#" << std::setw(width - 1) << "Vertex" << " ";
        for(; h_it!=h_end; ++h_it)
        {
          output_file << std::setw(width) << *h_it << " ";
        }
        output_file << std::endl;
      }
      for (std::tie(v_it, v_end) = boost::vertices(graph); v_it!=v_end; ++v_it)
      {
        output_file << std::setw(width) << graph[*v_it].name << " ";
        p_it = props.begin();
        for(; p_it!=p_end; ++p_it)
        {
          output_file << std::setw(width) << p_it->at(*v_it) << " ";
        }
        output_file << std::endl;
      }
    }
    else if(vertexID == notBGL::vertexID_num)
    {
      if(header.size() > 0)
      {
        output_file << "#" << std::setw(width - 1) << "Vertex" << " ";
        for(; h_it!=h_end; ++h_it)
        {
          output_file << std::setw(width) << *h_it << " ";
        }
        output_file << std::endl;
      }
      for (std::tie(v_it, v_end) = boost::vertices(graph); v_it!=v_end; ++v_it)
      {
        output_file << std::setw(width) << graph[*v_it].num << " ";
        p_it = props.begin();
        for(; p_it!=p_end; ++p_it)
        {
          output_file << std::setw(width) << p_it->at(*v_it) << " ";
        }
        output_file << std::endl;
      }
    }
    else if(vertexID == notBGL::vertexID_none)
    {
      if(header.size() > 0)
      {
        output_file << "#" << std::setw(width - 1) << *h_it << " ";
        for(++h_it; h_it!=h_end; ++h_it)
        {
          output_file << std::setw(width) << *h_it << " ";
        }
        output_file << std::endl;
      }
      for (std::tie(v_it, v_end) = boost::vertices(graph); v_it!=v_end; ++v_it)
      {
        p_it = props.begin();
        for(; p_it!=p_end; ++p_it)
        {
          output_file << std::setw(width) << p_it->at(*v_it) << " ";
        }
        output_file << std::endl;
      }
    }
    else
    {
      std::cerr << "Unknown vertex identifier type." << std::endl;
      std::terminate();
    }
  }
  // Closes the stream.
  output_file.close();
}


// ================================================================================================
// ================================================================================================
template<typename graph_t, typename vertex_t>
void notBGL::save_vertices_properties(std::string filename, graph_t &graph, std::initializer_list<std::map<vertex_t, double> > props, std::initializer_list<std::string> header, vertex_identifier_t vertexID, unsigned int width)
{
  notBGL::save_vertices_properties(filename, graph, props, vertexID, width, header);
}


// ================================================================================================
// ================================================================================================
template<typename graph_t, typename edge_t>
void notBGL::save_edges_properties(std::string filename, graph_t &graph, std::initializer_list<std::map<edge_t, double> > props, vertex_identifier_t vertexID, unsigned int width, std::initializer_list<std::string> header)
{
  // Stream objects.
  std::ofstream output_file;
  // Opens the stream and terminates if the operation did not succeed.
  output_file.open(filename.c_str(), std::ios_base::out);
  if( !output_file.is_open() )
  {
    std::cerr << "Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }
  else
  {
    // Iterators over the edges of the graph.
    typename boost::graph_traits<graph_t>::edge_iterator e_it, e_end;
    // Iterators over the std::string in "header".
    typename std::initializer_list<std::string>::iterator h_it(header.begin()), h_end(header.end());
    // Iterators over the std::map<edge_t, double> in "props".
    typename std::initializer_list<std::map<edge_t, double> >::iterator p_it, p_end(props.end());

    // Prints the vertex properties.
    if(vertexID == notBGL::vertexID_name)
    {
      if(header.size() > 0)
      {
        output_file << "#" << std::setw(width - 1) << "Vertex1" << " ";
        output_file << std::setw(width) << "Vertex2" << " ";
        for(; h_it!=h_end; ++h_it)
        {
          output_file << std::setw(width) << *h_it << " ";
        }
        output_file << std::endl;
      }
      for (std::tie(e_it, e_end) = boost::edges(graph); e_it!=e_end; ++e_it)
      {
        output_file << std::setw(width) << graph[boost::source(*e_it, graph)].name << " ";
        output_file << std::setw(width) << graph[boost::target(*e_it, graph)].name << " ";
        p_it = props.begin();
        for(; p_it!=p_end; ++p_it)
        {
          output_file << std::setw(width) << p_it->at(*e_it) << " ";
        }
        output_file << std::endl;
      }
    }
    else if(vertexID == notBGL::vertexID_num)
    {
      if(header.size() > 0)
      {
        output_file << "#" << std::setw(width - 1) << "Vertex1" << " ";
        output_file << std::setw(width) << "Vertex2" << " ";
        for(; h_it!=h_end; ++h_it)
        {
          output_file << std::setw(width) << *h_it << " ";
        }
        output_file << std::endl;
      }
      for (std::tie(e_it, e_end) = boost::edges(graph); e_it!=e_end; ++e_it)
      {
        output_file << std::setw(width) << graph[boost::source(*e_it, graph)].num << " ";
        output_file << std::setw(width) << graph[boost::target(*e_it, graph)].num << " ";
        p_it = props.begin();
        for(; p_it!=p_end; ++p_it)
        {
          output_file << std::setw(width) << p_it->at(*e_it) << " ";
        }
        output_file << std::endl;
      }
    }
    else if(vertexID == notBGL::vertexID_none)
    {
      if(header.size() > 0)
      {
        output_file << "#" << std::setw(width - 1) << *h_it << " ";
        for(++h_it; h_it!=h_end; ++h_it)
        {
          output_file << std::setw(width) << *h_it << " ";
        }
        output_file << std::endl;
      }
      for (std::tie(e_it, e_end) = boost::edges(graph); e_it!=e_end; ++e_it)
      {
        p_it = props.begin();
        for(; p_it!=p_end; ++p_it)
        {
          output_file << std::setw(width) << p_it->at(*e_it) << " ";
        }
        output_file << std::endl;
      }
    }
    else
    {
      std::cerr << "Unknown vertex identifier type." << std::endl;
      std::terminate();
    }
  }
  // Closes the stream.
  output_file.close();
}


// ================================================================================================
// ================================================================================================
template<typename graph_t, typename edge_t>
void notBGL::save_edges_properties(std::string filename, graph_t &graph, std::initializer_list<std::map<edge_t, double> > props, std::initializer_list<std::string> header, vertex_identifier_t vertexID, unsigned int width)
{
  notBGL::save_edges_properties(filename, graph, props, vertexID, width, header);
}










// ================================================================================================
// ================================================================================================
// *** Module: Topology
// ================================================================================================


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::degrees(graph_t &graph)
{
  // Iterators.
  typename boost::graph_traits<graph_t>::vertex_iterator v_it, v_end;
  // std::map containing the degrees.
  std::map<typename graph_t::vertex_descriptor, double> Vertex2Degree;
  // Fills the container.
  for (std::tie(v_it, v_end) = boost::vertices(graph); v_it!=v_end; ++v_it)
  {
    Vertex2Degree[*v_it] = boost::out_degree(*v_it, graph);
  }
  // Returns the degrees.
  return Vertex2Degree;
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::in_degrees(graph_t &graph)
{
  // Iterators.
  typename boost::graph_traits<graph_t>::vertex_iterator v_it, v_end;
  // std::map containing the degrees.
  std::map<typename graph_t::vertex_descriptor, double> Vertex2Degree;
  // Fills the container.
  for (std::tie(v_it, v_end) = boost::vertices(graph); v_it!=v_end; ++v_it)
  {
    Vertex2Degree[*v_it] = boost::in_degree(*v_it, graph);
  }
  // Returns the degrees.
  return Vertex2Degree;
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::out_degrees(graph_t &graph)
{
  // Returns the out-degrees.
  return notBGL::degrees(graph);
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::survey_triangles(graph_t &g)
{
  // Typedef.
  typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
  // typedef typename boost::graph_traits<graph_t>::edge_descriptor edge_t;
  // Vector objects.
  std::vector<vertex_t> intersection;
  std::vector<std::vector<vertex_t> > triangles;
  // std::vector<std::vector<vertex_t> > triangles;
  // Set objects.
  std::set<vertex_t> neighbours_v1, neighbours_v2;
  // Iterator objects.
  typename std::vector<vertex_t>::iterator it;
  typename boost::graph_traits<graph_t>::vertex_iterator it1, end1;
  typename boost::graph_traits<graph_t>::adjacency_iterator it2, end2, it3, end3;
  // Vertex descriptor objects.
  vertex_t v3;
  // Objects containing the degree of vertices.
  typename boost::graph_traits<graph_t>::degree_size_type d1, d2;
  // Computes the intersection for the in- and out- neighbourhoods of each node.
  for(std::tie(it1, end1) = boost::vertices(g); it1 != end1; ++it1)
  {
    // Degree of vertex v1.
    d1 = out_degree(*it1, g);
    // Performs the calculation only if d1>1.
    if( d1 > 1 )
    {
      // Builds an ordered list of the neighbourhood of v1
      neighbours_v1.clear();
      std::tie(it2, end2) = boost::adjacent_vertices(*it1, g);
      neighbours_v1.insert(it2, end2);
      // Loops over the neighbours of vertex v1.
      for(std::tie(it2, end2) = boost::adjacent_vertices(*it1, g); it2!=end2; ++it2)
      {
        // Identity and degree of vertex 2.
        d2 = out_degree(*it2, g);
        // Performs the calculation only if d1>1 and if v2>v1 (ensures that each triangle is counted once).
        if( *it1 < *it2 && d2 > 1 )
        {        
          // Builds an ordered list of the neighbourhood of v2
          neighbours_v2.clear();
          for(std::tie(it3, end3) = boost::adjacent_vertices(*it2, g); it3 != end3; ++it3)
          {
            if(*it2 < *it3) // Ensures that triangles will be counted only once.
            {
              neighbours_v2.insert(*it3);
            }
          }
          // Identifies the triangles.
          intersection.clear();
          intersection.resize(std::min(d1, neighbours_v2.size()));
          it = std::set_intersection(neighbours_v1.begin(), neighbours_v1.end(), neighbours_v2.begin(), neighbours_v2.end(), intersection.begin());
          intersection.resize(it-intersection.begin());
          // Loops over the common neighbours of vertices v1 and v2.
          for(unsigned int n(0), nn(intersection.size()); n<nn; ++n)
          {
            // Identity and degree of vertex 2.
            v3 = intersection[n];
            // Adds the triangle to the list.
            triangles.push_back({*it1, *it2, intersection[n]});
          }
        }
      }
    }
  }
  // Returns the vector containing the triangles.
  return triangles;
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::survey_directed_triangles(graph_t &g)
{
  // Typedef.
  typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
  // Vector objects.
  // std::vector<vertex_t> triangle(3);
  std::vector<vertex_t> intersection;
  // List objects.
  std::list<std::vector<vertex_t> > triangles;
  // Set objects.
  std::set<vertex_t> neighbours_v1, neighbours_v2;
  // Iterator objects.
  typename std::vector<vertex_t>::iterator it;
  typename std::set<vertex_t>::iterator it5, end5;
  typename boost::graph_traits<graph_t>::vertex_iterator it1, end1;
  typename boost::graph_traits<graph_t>::adjacency_iterator it2, end2;
  typename graph_t::inv_adjacency_iterator it4, end4;
  // Objects containing the degree of vertices.
  typename boost::graph_traits<graph_t>::degree_size_type d1, d2;
  // Computes the intersection for the in- and out- neighbourhoods of each node.
  for(std::tie(it1, end1) = boost::vertices(g); it1 != end1; ++it1)
  {
    // Degree of vertex v1.
    d1 = out_degree(*it1, g);
    d1 += in_degree(*it1, g);
    // Performs the calculation only if d1>1.
    if( d1 > 1 )
    {
      // Builds an ordered list of the neighbourhood of v1
      neighbours_v1.clear();
      for(std::tie(it2, end2) = boost::adjacent_vertices(*it1, g); it2!=end2; ++it2)
      {
        neighbours_v1.insert(*it2);
      }
      for(std::tie(it4, end4) = boost::inv_adjacent_vertices(*it1, g); it4!=end4; ++it4)
      {
        neighbours_v1.insert(*it4);
      }
      // Loops over the neighbours of vertex v1.
      it5 = neighbours_v1.begin();
      end5 = neighbours_v1.end();
      for(; it5 != end5; ++it5)
      {
        // Identity and degree of vertex 2.
        d2 = out_degree(*it5, g);
        d2 += in_degree(*it5, g);
        // Performs the calculation only if d1>1 and if v2>v1 (ensures that each triangle is counted once).
        if( *it1 < *it5 && d2 > 1 )
        {        
          // Builds an ordered list of the neighbourhood of v2
          neighbours_v2.clear();
          for(std::tie(it2, end2) = boost::adjacent_vertices(*it5, g); it2 != end2; ++it2)
          {
            if(*it5 < *it2) // Ensures that triangles will be counted only once.
            {
              neighbours_v2.insert(*it2);
            }
          }
          for(std::tie(it4, end4) = boost::inv_adjacent_vertices(*it5, g); it4 != end4; ++it4)
          {
            if(*it5 < *it4) // Ensures that triangles will be counted only once.
            {
              neighbours_v2.insert(*it4);
            }
          }
          // Identifies the triangles.
          intersection.clear();
          intersection.resize(std::min(neighbours_v1.size(), neighbours_v2.size()));
          it = std::set_intersection(neighbours_v1.begin(), neighbours_v1.end(), neighbours_v2.begin(), neighbours_v2.end(), intersection.begin());
          intersection.resize(it-intersection.begin());
          // Loops over the common neighbours of vertices v1 and v2.
          // triangle[0] = *it1;
          // triangle[1] = *it5;
          for(unsigned int n(0), nn(intersection.size()); n<nn; ++n)
          {
            // triangle[2] = intersection[n];
            // Adds the triangle to the list.
            triangles.push_back({*it1, *it5, intersection[n]});
          }
        }
      }
    }
  }
  // Returns the vector containing the triangles.
  return triangles;
}


// ================================================================================================
// ================================================================================================
template<typename vector_t, typename graph_t>
auto notBGL::local_clustering_coefficients(vector_t &triangles, graph_t &graph)
{
  // Typedef.
  typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
  // Map objects.
  std::map<vertex_t, double> local_clustering_coefficient;
  // Iterator objects.
  typename std::map<vertex_t, double>::iterator it;
  // Counts the number of triangles each vertex is part of.
  bool inserted;
  vertex_t v;
  for(auto t : triangles)
  {
    for(unsigned int i(0); i<3; ++i)
    {
      v = t[i];
      std::tie(it, inserted) = local_clustering_coefficient.insert(std::make_pair(v, double(0)));
      if(inserted)
      {
        it->second = 1;
      }
      else
      {
        it->second += 1;
      }
    }
  }
  // Computes the local clustering coefficient.
  typename boost::graph_traits<graph_t>::vertex_iterator v_it, v_end;
  typename boost::graph_traits<graph_t>::degree_size_type d;
  for(std::tie(v_it, v_end) = boost::vertices(graph); v_it!=v_end; ++v_it)
  {
    d = boost::out_degree(*v_it, graph);
    if(d > 1)
    {
      local_clustering_coefficient[*v_it] /= d * (d - 1) / 2;
    }
    else
    {
      local_clustering_coefficient[*v_it] = 0;
    }
  }
  // Returns the std::map mapping a vertex descriptior to a value of local clustering.
  return local_clustering_coefficient;
}


// // ================================================================================================
// // ================================================================================================
template<typename map_t>
double notBGL::average_local_clustering_coefficient(map_t &local_clustering_coefficients)
{
  return notBGL::utilities::average_of_map(local_clustering_coefficients);
}


// // ================================================================================================
// // ================================================================================================
template<typename vector_t, typename graph_t>
double notBGL::global_clustering_coefficient(vector_t &triangles, graph_t &graph)
{
  // Iterators.
  typename boost::graph_traits<graph_t>::vertex_iterator v_it, v_end;
  // Number of triplets in the graph.
  double nb_triplets = 0;
  unsigned int degree;
  for (std::tie(v_it, v_end) = boost::vertices(graph); v_it!=v_end; ++v_it)
  {
    degree = boost::out_degree(*v_it, graph);
    nb_triplets += degree * (degree - 1) / 2;
  }
  // Returns the global coefficient of clustering.
  return 3 * triangles.size() / nb_triplets;
}


// ================================================================================================
// ================================================================================================
template<typename vector_t, typename graph_t>
auto notBGL::multiplicity(vector_t &triangles, graph_t &graph)
{
  // Iterators.
  typename vector_t::iterator it, end;
  it = triangles.begin();
  end = triangles.end();
  // Edge descriptor.
  typename graph_t::edge_descriptor e;
  // Variable.
  bool exists;
  // Initializes the std::map containing the multiplicities.
  std::map<typename graph_t::edge_descriptor, double> Edge2Multiplicity;
  typename boost::graph_traits<graph_t>::edge_iterator e_it, e_end;
  for(std::tie(e_it, e_end) = boost::edges(graph); e_it!=e_end; ++e_it)
  {
    Edge2Multiplicity[*e_it] = 0;
  }
  // Loops over the triangles.
  for(; it!=end; ++it)
  {
    // Gets the vertices in the triangle.
    for(unsigned int i(0); i<3; ++i)
    {
      for(unsigned int j(i+1); j<3; ++j)
      {
        // Gets the edge descriptor.
        std::tie(e, exists) = boost::edge(it->at(i), it->at(j), graph);
        if(exists)
        {
          Edge2Multiplicity[e] += 1;
        }
        else
        {
          std::cerr << "Edge in triangle does not exists." << std::endl;
          std::terminate();
        }
      }
    }
  }
  // Returns the multiplicities.
  return Edge2Multiplicity;
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::reciprocical_edge_pairs(graph_t &graph)
{
  // Vertex descriptor.
  typedef typename graph_t::vertex_descriptor vertex_t;
  // Vector objects.
  typename std::vector<vertex_t> intersection;
  // Set object.
  typename std::set<vertex_t> in_neighbours, out_neighbours;
  // Map object.
  typename std::map<vertex_t, double> Vertex2NbReciprocalEdges;
  // Iterator objects.
  typename std::vector<vertex_t>::iterator it;
  typename graph_t::vertex_iterator v_it, v_end;
  typename graph_t::inv_adjacency_iterator in_it, in_end;
  typename graph_t::adjacency_iterator out_it, out_end;
  // In and out degrees.
  typename graph_t::degree_size_type in_degree, out_degree;
  // // Computes the intersection for the in- and out- neighbourhoods of each node.
  for(std::tie(v_it, v_end) = boost::vertices(graph); v_it!=v_end; ++v_it)
  {
    // In-degree.
    in_degree = boost::in_degree(*v_it, graph);
    in_neighbours.clear();
    for(std::tie(in_it, in_end) = boost::inv_adjacent_vertices(*v_it, graph); in_it!=in_end; ++in_it)
    {
      in_neighbours.insert(*in_it);
    }
    // Out-degree.
    out_degree = boost::out_degree(*v_it, graph);
    out_neighbours.clear();
    for(std::tie(out_it, out_end) = boost::adjacent_vertices(*v_it, graph); out_it!=out_end; ++out_it)
    {
      out_neighbours.insert(*out_it);
    }
    // Counts the number of edges that are reciprocal.
    intersection.clear();
    intersection.resize(in_degree, out_degree);
    it = std::set_intersection(in_neighbours.begin(), in_neighbours.end(), out_neighbours.begin(), out_neighbours.end(), intersection.begin());
    intersection.resize(it-intersection.begin());
    Vertex2NbReciprocalEdges[*v_it] = intersection.size();
  }
  // Returns the std::map mapping the vertex descriptor to the number of reciprocical edges.
  return Vertex2NbReciprocalEdges;
}


// ================================================================================================
// ================================================================================================
  template<typename map_t, typename graph_t>
  auto notBGL::reciprocity(map_t &Vector2ReciprocicalPairs, graph_t &graph)
  {
    // Gets the number of vertices.
    double nb_vertices = boost::num_vertices(graph);
    // Gets the number of directed edges.
    double nb_edges = boost::num_edges(graph);
    // Gets the number of edges that have a twin in the other direction.
    double nb_reciprocal_edge_pairs = notBGL::utilities::sum_of_map(Vector2ReciprocicalPairs);
    // Computes the fraction of edges that are reciprocical (traditional definition).
    double r = nb_reciprocal_edge_pairs / nb_edges;
    // Computes the reciprocity as defined in doi:10.1103/PhysRevLett.93.268701.
    double a_bar = nb_edges / (nb_vertices * (nb_vertices - 1));
    double rho = (r - a_bar) / (1 - a_bar);
    // Returns the two definitions.
    return std::make_tuple(r, rho);
  }


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::betweenness_centrality(graph_t &graph)
{
  // This code is greatly inspired from the examples given in
  // - http://programmingexamples.net/wiki/CPP/Boost/BGL/BetweennessCentralityClustering
  // - http://liuweipingblog.cn/cpp/an-example-of-boost-betweenness-centrality/
  // - http://stackoverflow.com/questions/23260793/to-convert-internal-properties-in-boost
  //                                                      -graph-to-external-properties-container-i


  // *** Creates a property map that accumulates the betweenness centrality of each vertex (could
  //   std::map be used?).

  // Typedef of a vertex.
  typedef typename graph_t::vertex_descriptor vertex_t;
  // Typedef of an index map for the vertices.
  typedef typename boost::property_map< graph_t, boost::vertex_index_t>::type VertexIndexMap;
  // Gets the actual graph's vertex indexmap.
  VertexIndexMap vertexIndexMap = get(boost::vertex_index, graph);
  // Creates a container with size equal to number of vertices in "graph".
  std::vector<double> v_centrality_vec(num_vertices(graph), 0.0);
  // Creates the property map.
  boost::iterator_property_map< std::vector<double>::iterator, VertexIndexMap >
    v_centrality_map(v_centrality_vec.begin(), vertexIndexMap);


  // *** Creates a property map that accumulates the betweenness centrality of each edge. Since we
  //   typically do not use boost::vecS as the container for the edges, we need to define a
  //   numerical ID for each edge.

  // Typedef of an edge.
  typedef typename graph_t::edge_descriptor edge_t;
   // std::map used for convenient initialization
  typedef typename std::map<edge_t, int> StdEdgeIndexMap;
  StdEdgeIndexMap my_e_index;
  // associative property map needed for iterator property map-wrapper
  typedef boost::associative_property_map< StdEdgeIndexMap > EdgeIndexMap;
  EdgeIndexMap e_index(my_e_index);
  // We use setS as edge-container -> no automatic indices
  // -> Create and set it explicitly
  int i = 0;
  typename boost::graph_traits<graph_t>::edge_iterator e_it, e_end;
  for(std::tie(e_it, e_end) = boost::edges(graph); e_it != e_end; ++e_it)
  {
    my_e_index.insert(std::pair< edge_t, int >(*e_it, i));
    ++i;
  }
  // Define EdgeCentralityMap
  std::vector< double > e_centrality_vec(boost::num_edges(graph), 0.0);
  // Create the external property map
  boost::iterator_property_map< std::vector< double >::iterator, EdgeIndexMap >
    e_centrality_map(e_centrality_vec.begin(), e_index);
 

  // *** Calculates the vertex and edge centralites.
  boost::brandes_betweenness_centrality(graph, v_centrality_map, e_centrality_map);


  // *** Converts the property maps into std::map (may not be necessary).
  std::map<vertex_t, double> Vertex2BC;
  typename boost::graph_traits<graph_t>::vertex_iterator v_it, v_end;
  for(std::tie(v_it, v_end) = boost::vertices(graph); v_it != v_end; ++v_it)
  {
    Vertex2BC[*v_it] = v_centrality_map[*v_it];
  }
  std::map<edge_t, double> Edge2BC;
  // typename boost::graph_traits<graph_t>::edge_iterator e_it, e_end;
  for(std::tie(e_it, e_end) = boost::edges(graph); e_it != e_end; ++e_it)
  {
    Edge2BC[*e_it] = e_centrality_map[*e_it];
  }


  // *** Returns a tuple containing both Vertex2BC and Edge2BC.
  return std::make_tuple(Vertex2BC, Edge2BC);
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::connected_components(graph_t &graph)
{
  // Typedef of a vertex.
  typedef typename graph_t::vertex_descriptor vertex_t;
  // Maps and associative maps to store the component memberships.
  typedef std::map<vertex_t, int> map_t;
  map_t component;
  boost::associative_property_map<map_t> component_ass_map(component);
  // Identifies the components.
  unsigned int num_components = boost::connected_components(graph, component_ass_map);
  // Returns the number of components and the memberships.
  return std::make_tuple(component, num_components);
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::kcore_decomposition(graph_t &graph)
{
  // Typedef of a vertex.
  typedef typename graph_t::vertex_descriptor vertex_t;
  // std::map object containing the coreness of each vertex.
  typename std::map<vertex_t, double> Vertex2kCore;
  // std::vector object containing the effective degree of each vertex.
  std::vector<unsigned int> DegreeVec(boost::num_vertices(graph));
  // Builds the list of the degree of the vertices.
  typename graph_t::vertex_iterator v_it, v_end;
  typename graph_t::adjacency_iterator a_it, a_end;
  std::set<std::pair<unsigned int, vertex_t> > DegreeSet;
  for (std::tie(v_it, v_end) = boost::vertices(graph); v_it!=v_end; ++v_it)
  {
    DegreeSet.insert(std::make_pair(boost::out_degree(*v_it, graph), *v_it));
    DegreeVec[graph[*v_it].num] = boost::out_degree(*v_it, graph);
  }
  // Sets the coreness of each vertex using the algorithm by  Batagelj and Zaversnik.
  vertex_t v1, v2;
  unsigned int d1, d2, n2;
  typename std::set<std::pair<unsigned int, vertex_t> >::iterator m_it;
  while(DegreeSet.size() > 0)
  {
    // Sets the coreness of the first vertex of the list and removes it from the list.
    d1 = DegreeSet.begin()->first;
    v1 = DegreeSet.begin()->second;
    Vertex2kCore[v1] = d1;
    DegreeSet.erase(DegreeSet.begin());
    // Reduces the "effective" degree of its neighbours.
    for(std::tie(a_it, a_end) = boost::adjacent_vertices(v1, graph); a_it!=a_end; ++a_it)
    {
      n2 = graph[*a_it].num;
      d2 = DegreeVec[n2];
      m_it = DegreeSet.find(std::make_pair(d2, *a_it));
      if(m_it != DegreeSet.end())
      {
        // std::cout << graph[m_it->second].name << std::endl;
        d2 = m_it->first;
        if(d2 > d1)
        {
          v2 = m_it->second;
          DegreeVec[n2] = d2 - 1;
          DegreeSet.erase(m_it);
          DegreeSet.insert(std::make_pair(d2 - 1, v2));
        }
      }
    }
  }
  // Returns the coreness of the vertices.
  return Vertex2kCore;
}


// ================================================================================================
// ================================================================================================
template<typename vector_t, typename graph_t>
std::vector<double> notBGL::triangle_spectrum(vector_t &triangles, graph_t &graph)
{
  // Typedef.
  typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
  // Vector objects.
  std::vector<double> triangle_spectrum(7, double(0));
  // Defines the homomorphism groups.
  boost::multi_array<unsigned int, 3> config(boost::extents[3][3][3]);
  config[0][0][0] = 0;  config[1][0][0] = 1;  config[2][0][0] = 3;
  config[0][0][1] = 1;  config[1][0][1] = 1;  config[2][0][1] = 2;
  config[0][0][2] = 3;  config[1][0][2] = 4;  config[2][0][2] = 5;
  config[0][1][0] = 1;  config[1][1][0] = 1;  config[2][1][0] = 4;
  config[0][1][1] = 1;  config[1][1][1] = 0;  config[2][1][1] = 3;
  config[0][1][2] = 2;  config[1][1][2] = 3;  config[2][1][2] = 5;
  config[0][2][0] = 3;  config[1][2][0] = 2;  config[2][2][0] = 5;
  config[0][2][1] = 4;  config[1][2][1] = 3;  config[2][2][1] = 5;
  config[0][2][2] = 5;  config[1][2][2] = 5;  config[2][2][2] = 6;
  // Variables used to identify the homomorphism groups.
  unsigned int edge01, edge12, edge20;
  bool is_edge1, is_edge2;
  typename graph_t::edge_descriptor e;
  // Classifies the triangles.
  for(auto t : triangles)
  {
    // Edge 0 <--> 1.
    std::tie(e, is_edge1) = boost::edge(t[0], t[1], graph);
    std::tie(e, is_edge2) = boost::edge(t[1], t[0], graph);
    if(is_edge1)
    {
      if(is_edge2)
        edge01 = 2;
      else
        edge01 = 0;
    }
    else
      edge01 = 1;
    // Edge 1 <--> 2.
    std::tie(e, is_edge1) = boost::edge(t[1], t[2], graph);
    std::tie(e, is_edge2) = boost::edge(t[2], t[1], graph);
    if(is_edge1)
    {
      if(is_edge2)
        edge12 = 2;
      else
        edge12 = 0;
    }
    else
      edge12 = 1;
    // Edge 2 <--> 0.
    std::tie(e, is_edge1) = boost::edge(t[2], t[0], graph);
    std::tie(e, is_edge2) = boost::edge(t[0], t[2], graph);
    if(is_edge1)
    {
      if(is_edge2)
        edge20 = 2;
      else
        edge20 = 0;
    }
    else
      edge20 = 1;

    // Compiles the triangle class.
    triangle_spectrum[config[edge01][edge12][edge20]] += 1;

  }
  // Returns the std::map mapping a vertex descriptior to a value of local clustering.
  return triangle_spectrum;
}










// ================================================================================================
// ================================================================================================
// *** Module: Weighted organization
// ================================================================================================


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::strengths(graph_t &graph)
{
  // std::map containing the weights.
  std::map<typename graph_t::vertex_descriptor, double> Vertex2Strength;
  // Iterators over the edges of the graph.
  typename boost::graph_traits<graph_t>::vertex_iterator v_it, v_end;
  // Copies the weights.
  for(std::tie(v_it, v_end) = boost::vertices(graph); v_it!=v_end; ++v_it)
  {
    Vertex2Strength[*v_it] = graph[*v_it].strength;
  }
  // Returns the weights.
  return Vertex2Strength;
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::weights(graph_t &graph)
{
  // Internal property map for the weights.
  typename boost::property_map<graph_t, boost::edge_weight_t>::type weight_map = get(boost::edge_weight, graph);
  // std::map containing the weights.
  std::map<typename graph_t::edge_descriptor, double> Edge2Weight;
  // Iterators over the edges of the graph.
  typename boost::graph_traits<graph_t>::edge_iterator e_it, e_end;
  // Copies the weights.
  for(std::tie(e_it, e_end) = boost::edges(graph); e_it!=e_end; ++e_it)
  {
    Edge2Weight[*e_it] = get(weight_map, *e_it);
  }
  // Returns the weights.
  return Edge2Weight;
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::disparity(graph_t &graph)
{
  // Typedef.
  typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
  // Objects related to the weights.
  typename boost::property_map<graph_t, boost::edge_weight_t>::type weight_map = get(boost::edge_weight, graph);
  // Iterators.
  typename boost::graph_traits<graph_t>::vertex_iterator v_it, v_end;
  typename boost::graph_traits<graph_t>::out_edge_iterator e_it, e_end;
  // Map objects.
  std::map<vertex_t, double> disparity;
  // Computes the disparity for each vertex.
  bool inserted;
  double tmp;
  for(std::tie(v_it, v_end) = boost::vertices(graph); v_it!=v_end; ++v_it)
  {
    // Initializes the value of the disparity.
     tmp = 0;
    // Runs over every neighbors.
    for(std::tie(e_it, e_end) = boost::out_edges(*v_it, graph); e_it!=e_end; ++e_it)
    {
      tmp += std::pow(get(weight_map, *e_it), 2);
      // std::cout << get(weight_map, *e_it) << std::endl;
    }
    if(graph[*v_it].strength > 0)
    {
      disparity[*v_it] = tmp / std::pow(graph[*v_it].strength, 2);
    }
  }
  // Returns the disparity.
  return disparity;
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::load_weighted_edgelist(std::string filename, graph_t &graph)
{
  // Vertex descriptors.
  typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
  vertex_t v1, v2;
  // Objects related to edges.
  typedef typename graph_t::edge_property_type weight_t;
  double weight;
  // Stream objects.
  std::ifstream edgelist_file;
  std::stringstream one_line;
  // String objects.
  std::string full_line, name1_str, name2_str, weight_str;
  // Integer objects.
  unsigned int node_cnt(0);
  // Map objects to assure the uniqueness of vertices.
  std::map< std::string, vertex_t > Name2Vertex;
  // Iterator objects.
  typename std::map< std::string, vertex_t >::iterator name_it;
  // Opens the stream and terminates if the operation did not succeed.
  edgelist_file.open(filename.c_str(), std::ios_base::in);
  if( !edgelist_file.is_open() )
  {
    std::cerr << "Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }
  else
  {
    // Reads the edgelist and registers the nodes and the edges.
    while( !edgelist_file.eof() )
    {
      // Reads a line of the edgelist.
      std::getline(edgelist_file,full_line); edgelist_file >> std::ws;
      one_line.str(full_line); one_line >> std::ws;
      one_line >> name1_str >> std::ws;
      one_line >> name2_str >> std::ws;
      one_line >> weight_str >> std::ws;
      one_line.clear();
      // Is name1 new?
      name_it = Name2Vertex.find(name1_str);
      if( name_it == Name2Vertex.end() )
      {
        v1 = boost::add_vertex(graph);
        Name2Vertex[name1_str] = v1;
        graph[v1].name = name1_str;
        graph[v1].num = node_cnt;
        ++node_cnt;
      }
      else
      {
        v1 = name_it->second;
      }
      // Is name2 new?
      name_it = Name2Vertex.find(name2_str);
      if( name_it == Name2Vertex.end() )
      {
        v2 = boost::add_vertex(graph);
        Name2Vertex[name2_str] = v2;
        graph[v2].name = name2_str;
        graph[v2].num = node_cnt;
        ++node_cnt;
      }
      else
      {
        v2 = name_it->second;
      }
      // Creates the edge if it does not exist (forbids multiple edges and self-loops).
      if(v1 != v2) // not a self-loop.
      {
        if(boost::edge(v1,v2,graph).second == false) // not a multiple edge
        {
          weight = std::stod(weight_str);
          graph[v1].strength += weight;
          graph[v2].strength += weight;
          boost::add_edge(v1, v2, weight_t(weight), graph);
        }
      }
    }
  }
  // Closes the stream.
  edgelist_file.close();
  // Returns the maps "Name2Vertex" to help add other properties to vertices.
  return Name2Vertex;
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
void notBGL::save_weighted_edgelist(std::string filename, graph_t &graph, vertex_identifier_t vertexID = vertexID_name, unsigned int width = 10)
{
  // Extracts the weights of the edges.
  auto Edge2Weight = notBGL::weights(graph);
  // Writes the weights into a file.
  notBGL::save_edges_properties(filename, graph, {Edge2Weight}, vertexID, width);
}










// ================================================================================================
// ================================================================================================
// *** Module: Geometry
// ================================================================================================


// ================================================================================================
// ================================================================================================
template<typename graph_t, typename map_t>
auto notBGL::load_coordinates(std::string filename, graph_t &graph, map_t &Name2Vertex)
{
  // Vertex descriptors.
  typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
  vertex_t v1;
  // Stream objects.
  std::ifstream positions_file;
  std::stringstream one_line;
  // String objects.
  std::string full_line, name1_str, pos1_str;
  // Map objects to assure the uniqueness of vertices.
  std::map< vertex_t, std::vector<double> > Vertex2Position;
  // Iterator objects.
  typename map_t::iterator name_it;
  // Vector objects.
  std::vector<double> position;
  // Variable.
  unsigned int N = 0;
  // Opens the stream and terminates if the operation did not succeed.
  positions_file.open(filename.c_str(), std::ios_base::in);
  if( !positions_file.is_open() )
  {
    std::cerr << "Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }
  else
  {
    // Reads the first line of the file to count the number of dimensions.
    std::getline(positions_file,full_line);
    one_line.str(full_line); one_line >> std::ws;
    one_line >> name1_str >> std::ws;
    // Skips a line of comment.
    while(name1_str == "#")
    {
      one_line.clear();
      std::getline(positions_file,full_line);
      one_line.str(full_line); one_line >> std::ws;
      one_line >> name1_str >> std::ws;
    }
    while(!one_line.eof())
    {
      one_line >> pos1_str >> std::ws;
      ++N;
    }
    one_line.clear();
    position.resize(N);
    // Resets the file stream to the initial position.
    positions_file.seekg(0, std::ios::beg);
    // Reads each line of the positions file.
    while( !positions_file.eof() )
    {
      // Reads a line of the positions file.
      std::getline(positions_file,full_line); positions_file >> std::ws;
      one_line.str(full_line); one_line >> std::ws;
      one_line >> name1_str >> std::ws;
      // Skips a line of comment.
      if(name1_str == "#")
      {
        one_line.clear();
        continue;
      }
      for(unsigned int i(0); i<N; ++i)
      {
        one_line >> pos1_str >> std::ws;
        position[i] = std::stod(pos1_str);
      }
      one_line.clear();
      // Gets the iterator corresponding to the vertex in the Name2Vertex map.
      name_it = Name2Vertex.find(name1_str);
      // If the vertex is not present in the edgelist (has a degree equal to 0)
      if( name_it == Name2Vertex.end() )
      {
        vertex_t v1 = boost::add_vertex(graph);
        Name2Vertex[name1_str] = v1;
        graph[v1].name = name1_str;
        graph[v1].num = boost::num_vertices(graph);
        Vertex2Position[v1] = position;
      }
      // Register the position of the vertex if it is present in the edgelist.
      else
      {
        Vertex2Position[name_it->second] = position;
      }
    }
  }
  // Closes the stream.
  positions_file.close();
  // Returns the maps "Vertex2Position".
  return Vertex2Position;
}

// ================================================================================================
// ================================================================================================
double notBGL::hyperbolic_distance(std::vector<double> &x1, std::vector<double> x2, double zeta)
{
  // Uses the hyperbolic law of cosines to compute the distance between to points in the
  // hyperbolic disk.
  double tmp = 3.14159265358979 - std::fabs(3.14159265358979 - std::fabs(x1[1] - x2[1]));
// if(tmp > 3.14159265358979)
//   std::cout << tmp << std::endl;
  tmp = std::sinh(zeta * x1[0]) * std::sinh(zeta * x2[0]) * std::cos(tmp);
  tmp = std::cosh(zeta * x1[0]) * std::cosh(zeta * x2[0]) - tmp;
  if(tmp < 1)
    tmp = 1;
  return std::acosh(tmp) / zeta;
}

// ================================================================================================
// ================================================================================================
double notBGL::euclidean_distance(std::vector<double> &x1, std::vector<double> x2)
{
  // Subtracts the two vectors.
  std::transform(x1.begin(), x1.end(), x2.begin(), x2.begin(), std::minus<double>());
  // Returns the square root of the inner product of the two vectors.
  return std::sqrt( std::inner_product(x2.begin(), x2.end(), x2.begin(), 0.0) );
}

#endif // NOTBGL_HPP_INCLUDED
