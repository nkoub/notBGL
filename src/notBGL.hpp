/**
 * @file notBGL.hpp
 *
 * @brief      Source code of the notBGL library.
 *
 *             This file contains the complete source code of the notBGL
 *             library. While having one single file is not the most clear and
 *             organized choice for source codes, it has been
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
#include <map>               // std::map
#include <numeric>           // std::inner_product
#include <queue>             // std::queue
#include <sstream>           // std::stringstream
#include <string>            // std::string, std::stod
#include <tuple>             // std::tie, std::make_tuple
#include <vector>            // std::vector
// Boost Graph Library
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>

// All the methods are under the namespace "notBGL"
namespace notBGL
{

  #ifndef DOXYGEN_EXCLUDED
    // Minimal vertex property structure used in the generic undirected unweighted graph type (not documented).
    struct MinimalVertexProp  { std::string name; unsigned int id; };
    struct MinimalWeightedVertexProp  { std::string name; unsigned int id; double strength = 0; };
  #endif // DOXYGEN_EXCLUDED

  /// Generic undirected unweighted graph type.
  typedef boost::adjacency_list<boost::listS,
                                boost::vecS,
                                boost::undirectedS,
                                MinimalVertexProp,
                                boost::no_property > UndirectedGraph_t;

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
   *                       problem.
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
   *
   * @tparam     graph_t   boost::adjacency_list
   * @tparam     vertex_t  boost::vertex_descriptor
   *
   * @see        save_edges_properties()
   *
   * @ingroup    IO
   */
  template<typename graph_t, typename vertex_t>
  void save_vertices_properties(std::string filename, graph_t &graph, std::initializer_list<std::map<vertex_t, double> > props, vertex_identifier_t vertexID = vertexID_name, unsigned int width = 10);

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
   *
   * @tparam     graph_t   boost::adjacency_list
   * @tparam     vertex_t  boost::edge_descriptor
   *
   * @see        save_vertices_properties()
   *
   * @ingroup    IO
   */
  template<typename graph_t, typename edge_t>
  void save_edges_properties(std::string filename, graph_t &graph, std::initializer_list<std::map<edge_t, double> > props, vertex_identifier_t vertexID = vertexID_name, unsigned int width = 10);





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
  template <typename vector_t, typename graph_t>
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
      double average_of_map(map_t &map);
    }
  #endif // DOXYGEN_EXCLUDED

}

// // ================================================================================================
// // ================================================================================================
#ifndef DOXYGEN_EXCLUDED
template<typename map_t>
double notBGL::utilities::average_of_map(map_t &map)
{
  double sum = 0;
  for(auto el : map)
  {
    sum += el.second;
  }
  return sum / map.size();
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
      one_line >> name2_str >> std::ws;
      one_line.clear();
      // Is name1 new?
      name_it = Name2Vertex.find(name1_str);
      if( name_it == Name2Vertex.end() )
      {
        v1 = boost::add_vertex(g);
        Name2Vertex[name1_str] = v1;
        g[v1].name = name1_str;
        g[v1].id = node_cnt;
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
        g[v2].id = node_cnt;
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
        edgelist_file << g[source(*it, g)].id << " "
                      << g[target(*it, g)].id << std::endl;
      }
    }
  }
  // Closes the stream.
  edgelist_file.close();
}


// ================================================================================================
// ================================================================================================
template<typename graph_t, typename vertex_t>
void notBGL::save_vertices_properties(std::string filename, graph_t &graph, std::initializer_list<std::map<vertex_t, double> > props, vertex_identifier_t vertexID = vertexID_name, unsigned int width = 10)
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
    // Iterators over the std::map<vertex_t, double> in "props".
    typename std::initializer_list<std::map<vertex_t, double> >::iterator p_it, p_end(props.end());

    // Prints the vertex properties.
    if(vertexID == notBGL::vertexID_name)
    {
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
      for (std::tie(v_it, v_end) = boost::vertices(graph); v_it!=v_end; ++v_it)
      {
        output_file << std::setw(width) << graph[*v_it].id << " ";
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
template<typename graph_t, typename edge_t>
void notBGL::save_edges_properties(std::string filename, graph_t &graph, std::initializer_list<std::map<edge_t, double> > props, vertex_identifier_t vertexID = vertexID_name, unsigned int width = 10)
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
    // Iterators over the std::map<edge_t, double> in "props".
    typename std::initializer_list<std::map<edge_t, double> >::iterator p_it, p_end(props.end());

    // Prints the vertex properties.
    if(vertexID == notBGL::vertexID_name)
    {
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
      for (std::tie(e_it, e_end) = boost::edges(graph); e_it!=e_end; ++e_it)
      {
        output_file << std::setw(width) << graph[boost::source(*e_it, graph)].id << " ";
        output_file << std::setw(width) << graph[boost::target(*e_it, graph)].id << " ";
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
  for(std::tie(it1, end1) = vertices(g); it1 != end1; ++it1)
  {
    // Degree of vertex v1.
    d1 = out_degree(*it1, g);
    // Performs the calculation only if d1>1.
    if( d1 > 1 )
    {
      // Builds an ordered list of the neighbourhood of v1
      neighbours_v1.clear();
      std::tie(it2, end2) = adjacent_vertices(*it1, g);
      neighbours_v1.insert(it2, end2);
      // Loops over the neighbours of vertex v1.
      for(std::tie(it2, end2) = adjacent_vertices(*it1, g); it2!=end2; ++it2)
      {
        // Identity and degree of vertex 2.
        d2 = out_degree(*it2, g);
        // Performs the calculation only if d1>1 and if v2>v1 (ensures that each triangle is counted once).
        if( *it1 < *it2 && d2 > 1 )
        {        
          // Builds an ordered list of the neighbourhood of v2
          neighbours_v2.clear();
          for(std::tie(it3, end3) = adjacent_vertices(*it2, g); it3 != end3; ++it3)
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
        graph[v1].id = node_cnt;
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
        graph[v2].id = node_cnt;
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
  // Vectro objects.
  std::vector<double> position();
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
    while(!one_line.eof())
    {
      ++N;
    }
    // Resets the file stream to the initial position.
    positions_file.seekg(0, std::ios::beg);
    // Reads each line of the positions file.
    while( !positions_file.eof() )
    {
      // Reads a line of the positions file.
      std::getline(positions_file,full_line); positions_file >> std::ws;
      one_line.str(full_line); one_line >> std::ws;
      one_line >> name1_str >> std::ws;
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
        graph[v1].id = boost::num_vertices(graph);
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
if(tmp > 3.14159265358979)
  std::cout << tmp << std::endl;
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
