/**
* @brief  Sets of methods to facilitate the manipulation of Boost Graph Library objects.
* @author Antoine Allard (<a href="http://antoineallard.info">antoineallard.info</a>)
* @date   March 2016
* @bug    No known bugs.
* @todo   Complete the documentation using doxygen.
*/

#ifndef NOTBGL_HPP_INCLUDED
#define NOTBGL_HPP_INCLUDED

// Standard Template Library
#include <algorithm>     // std::transform
#include <cmath>         // std::sqrt, std::pow, std::cosh, std::sinh, std::acosh, std::fabs
#include <cstdlib>       // std::terminate
#include <fstream>       // std::ifstream, std::ofstream
#include <functional>    // std::minus
#include <iomanip>       // std::setw
#include <iostream>      // std::cerr
#include <limits>        // std::numeric_limits
#include <map>           // std::map
#include <numeric>       // std::inner_product
#include <queue>         // std::queue
#include <sstream>       // std::stringstream
#include <string>        // std::string, std::stod
#include <tuple>         // std::tie, std::make_tuple
#include <vector>        // std::vector
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
  auto load_weighted_edgelist(std::string filename, graph_t &g);

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
  void save_weighted_edgelist(std::string filename, graph_t &g, bool write_names = true);

  /**
   * @brief      { function_description }
   *
   * @param[in]  filename           The filename
   * @param      g                  { parameter_description }
   * @param[in]  vertex_identifier  The vertex identifier
   *
   * @tparam     graph_t            { description }
   *
   * @ingroup    IO
   *
   */
  template<typename graph_t>
  void save_degree_sequence(std::string filename, graph_t &g, std::string vertex_identifier = "None");

  /**
   * @brief      { function_description }
   *
   * @param[in]  filename           The filename
   * @param      g                  { parameter_description }
   * @param[in]  vertex_identifier  The vertex identifier
   *
   * @tparam     graph_t            { description }
   *
   * @ingroup    IO
   */
  template<typename graph_t>
  void save_strength_sequence(std::string filename, graph_t &g, std::string vertex_identifier = "None");

  /**
   * @brief      { function_description }
   *
   * @param[in]  filename           The filename
   * @param      clustering         The clustering
   * @param      g                  { parameter_description }
   * @param[in]  vertex_identifier  The vertex identifier
   *
   * @tparam     map_t              { description }
   * @tparam     graph_t            { description }
   *
   * @ingroup    IO
   */
  template<typename map_t, typename graph_t>
  void save_local_clustering_coefficents_sequence(std::string filename, map_t &clustering, graph_t &g, std::string vertex_identifier = "None");

  /**
   * @brief      { function_description }
   *
   * @param[in]  filename    The filename
   * @param      g           { parameter_description }
   * @param[in]  clustering  The clustering
   * @param[in]  disparity   The disparity
   *
   * @tparam     graph_t     { description }
   * @tparam     map_t       { description }
   *
   * @ingroup    IO
   */
  template<typename graph_t, typename map_t>
  void save_vertices_properties(std::string filename, graph_t &g, map_t clustering, map_t disparity);





  /** ---------------------------------------------------------------------------------------------
   * @defgroup   topo Topology
   * @brief      This group contains the functions related to the
   *             characterization of the topology of graphs (e.g., clustering,
   *             betweenness centrality).
   */

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
   * @brief      Get the local clustering coefficient.
   *
   * @param[in]  triangles  The triangles
   * @param      g          { parameter_description }
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
  auto local_clustering_coefficient(vector_t triangles, graph_t &g);

  /**
   * @brief      Computes the betweenness centrality of the vertices and the
   *             edges of a graph.
   *
   * @param      graph    The graph object.
   *
   * @tparam     graph_t  A boost::adjacency_list type.
   *
   * @return     { description_of_the_return_value }
   *
   * @ingroup    topo
   */
  template<typename graph_t>
  auto betweenness_centrality(graph_t &graph);
  /**
   * @example    betweenness_centrality.cpp
   * @brief      This is an example of how to use the betweenness_centrality()
   *             function. See comments in the code and the documentation for
   *             further details.
   */




  /** ---------------------------------------------------------------------------------------------
   * @defgroup   routing Routing
   * @brief      .
   */

  /**
   * @brief      { function_description }
   *
   * @param      v_source  The v source
   * @param      g         { parameter_description }
   *
   * @tparam     vertex_t  { description }
   * @tparam     graph_t   { description }
   *
   * @return     { description_of_the_return_value }
   *
   * @ingroup    routing
   */
  template<typename vertex_t, typename graph_t>
  auto shortest_path_lengths_from_source(vertex_t &v_source, graph_t &g);





  /** ---------------------------------------------------------------------------------------------
   * @defgroup   geo Geometry
   * @brief      This group contains the functions related to graphs embedded in
   *             a geometric space (e.g., loading coordinates, calculating
   *             distances, greedy routing).
   */

  /**
   * @brief      { function_description }
   *
   * @param[in]  filename     The filename
   * @param      g            { parameter_description }
   * @param      Name2Vertex  The name2 vertex
   * @param[in]  N            { parameter_description }
   *
   * @tparam     graph_t      { description }
   * @tparam     map_t        { description }
   *
   * @return     { description_of_the_return_value }
   *
   * @ingroup    geo
   */
  template<typename graph_t, typename map_t>
  auto load_coordinates(std::string filename, graph_t &g, map_t &Name2Vertex, const unsigned int N = 3);

  /**
   * @brief      { function_description }
   *
   * @param      x1    { parameter_description }
   * @param[in]  x2    { parameter_description }
   *
   * @return     { description_of_the_return_value }
   *
   * @ingroup    geo
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
   *
   * @ingroup    geo
   */
  double hyperbolic_distance(std::vector<double> &x1, std::vector<double> x2, double zeta);


  


  /** ---------------------------------------------------------------------------------------------
   * @defgroup   weights Weighted organization
   * @brief      This group contains the functions related to graphs embedded in
   *             a geometric space (e.g., loading coordinates, calculating
   *             distances, greedy routing).
   */

  /**
   * @brief      { function_description }
   *
   * @param      g        { parameter_description }
   *
   * @tparam     graph_t  { description }
   *
   * @return     { description_of_the_return_value }
   *
   * @ingroup    weights
   */
  template<typename graph_t>
  auto disparity(graph_t &g);

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
template<typename graph_t, typename map_t>
void notBGL::save_vertices_properties(std::string filename, graph_t &g, map_t clustering, map_t disparity)
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
    typename boost::graph_traits<graph_t>::vertex_iterator it, end;
    for (std::tie(it, end) = vertices(g); it!=end; ++it)
    {
      output_file << std::setw(5) << g[*it].name              << " "
                  << std::setw(7)<< boost::out_degree(*it, g) << " "
                  << std::setw(7)<< g[*it].strength           << " "
                  << std::setw(9)<< clustering[*it]           << " "
                  << std::setw(9)<< disparity[*it]            << " "
                  << std::endl;
    }
  }
  // Closes the stream.
  output_file.close();
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
auto notBGL::disparity(graph_t &g)
{
  // Typedef.
  typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
  // Objects related to the weights.
  typename boost::property_map<graph_t, boost::edge_weight_t>::type weight_map = get(boost::edge_weight, g);
  // Iterators.
  typename boost::graph_traits<graph_t>::vertex_iterator v_it, v_end;
  typename boost::graph_traits<graph_t>::out_edge_iterator e_it, e_end;
  // Map objects.
  std::map<vertex_t, double> disparity;
  // Computes the disparity for each vertex.
  bool inserted;
  double tmp;
  for(std::tie(v_it, v_end) = boost::vertices(g); v_it!=v_end; ++v_it)
  {
    // Initializes the value of the disparity.
     tmp = 0;
    // Runs over every neighbors.
    for(std::tie(e_it, e_end) = boost::out_edges(*v_it, g); e_it!=e_end; ++e_it)
    {
      tmp += std::pow(get(weight_map, *e_it), 2);
      // std::cout << get(weight_map, *e_it) << std::endl;
    }
    if(g[*v_it].strength > 0)
    {
      disparity[*v_it] = tmp / std::pow(g[*v_it].strength, 2);
    }
  }
  //
  return disparity;
}

// ================================================================================================
// ================================================================================================
template<typename vector_t, typename graph_t>
auto notBGL::local_clustering_coefficient(vector_t triangles, graph_t &g)
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
  for(std::tie(v_it, v_end) = boost::vertices(g); v_it!=v_end; ++v_it)
  {
    d = boost::out_degree(*v_it, g);
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
template<typename vertex_t, typename graph_t>
auto shortest_path_lengths_from_source(vertex_t &v_source, graph_t &g)
{
  // Vector objects.
  std::vector<unsigned int> topological_distance(boost::num_vertices(g), 0);
  std::vector<bool> visited(boost::num_vertices(g), false);
  // Set objects.
  std::queue<vertex_t> to_visit;
  // Variables.
  unsigned int v_id;
  unsigned int nb_hop;
  // 
  v_id = g[v_source].id;
  to_visit.push(v_source);
  visited[v_id] = true;
  topological_distance[v_id] = nb_hop;
  //
  typename boost::graph_traits<graph_t>::adjacency_iterator it, end;
  while(!to_visit.empty())
  {
    nb_hop = topological_distance[g[to_visit.front()].id] + 1;
    for(std::tie(it, end) = adjacent_vertices(to_visit.front(), g); it != end; ++it)
    {
      v_id = g[*it].id;
      if(!visited[v_id])
      {
        topological_distance[v_id] = nb_hop;
        visited[v_id] = true;
        to_visit.push(*it);
      }
    }
    to_visit.pop();
  }
  //
  return topological_distance;
}

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
auto notBGL::load_weighted_edgelist(std::string filename, graph_t &g)
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
          weight = std::stod(weight_str);
          g[v1].strength += weight;
          g[v2].strength += weight;
          boost::add_edge(v1, v2, weight_t(weight), g);
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
template<typename graph_t>
void notBGL::save_degree_sequence(std::string filename, graph_t &g, std::string vertex_identifier = "None")
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
    typename boost::graph_traits<graph_t>::vertex_iterator it, end;
    if(vertex_identifier == "None")
    {
      // Writes the degree sequence, one vertex per line.
      for (std::tie(it, end) = vertices(g); it!=end; ++it)
      {
        output_file << out_degree(*it, g) << std::endl;
      }
    }
    else if(vertex_identifier == "Name")
    {
      // Writes the degree sequence, one vertex per line.
      for (std::tie(it, end) = vertices(g); it!=end; ++it)
      {
        output_file << g[*it].name        << " "
                    << out_degree(*it, g) << std::endl;
      }
    }
    else if(vertex_identifier == "ID")
    {
      // Writes the degree sequence, one vertex per line.
      for (std::tie(it, end) = vertices(g); it!=end; ++it)
      {
        output_file << g[*it].id          << " "
                    << out_degree(*it, g) << std::endl;
      }
    }
    else
    {
      std::cerr << "Unknown vertex identifier option: " << vertex_identifier << "." << std::endl;
      std::terminate();
    }
  }
  // Closes the stream.
  output_file.close();
}

// ================================================================================================
// ================================================================================================
template<typename map_t, typename graph_t>
void notBGL::save_local_clustering_coefficents_sequence(std::string filename, map_t &clustering, graph_t &g, std::string vertex_identifier = "None")
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
    // Writes the sequence, one vertex per line.
    if(vertex_identifier == "None")
    {
      for(auto v : clustering)
      {
        output_file << v.second << std::endl;
      }
    }
    else if(vertex_identifier == "Name")
    {
      for(auto v : clustering)
      {
        output_file << g[v.first].name << " "
                    << v.second << std::endl;
      }
    }
    else if(vertex_identifier == "ID")
    {
      for(auto v : clustering)
      {
        output_file << g[v.first].id << " "
                    << v.second << std::endl;
      }
    }
    else
    {
      std::cerr << "Unknown vertex identifier option: " << vertex_identifier << "." << std::endl;
      std::terminate();
    }
  }
  // Closes the stream.
  output_file.close();
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
void notBGL::save_strength_sequence(std::string filename, graph_t &g, std::string vertex_identifier = "None")
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
    typename boost::graph_traits<graph_t>::vertex_iterator it, end;
    if(vertex_identifier == "None")
    {
      // Writes the degree sequence, one vertex per line.
      for (std::tie(it, end) = vertices(g); it!=end; ++it)
      {
        output_file << g[*it].strength << std::endl;
      }
    }
    else if(vertex_identifier == "Name")
    {
      // Writes the degree sequence, one vertex per line.
      for (std::tie(it, end) = vertices(g); it!=end; ++it)
      {
        output_file << g[*it].name     << " "
                    << g[*it].strength << std::endl;
      }
    }
    else if(vertex_identifier == "ID")
    {
      // Writes the degree sequence, one vertex per line.
      for (std::tie(it, end) = vertices(g); it!=end; ++it)
      {
        output_file << g[*it].id       << " "
                    << g[*it].strength << std::endl;
      }
    }
    else
    {
      std::cerr << "Unknown vertex identifier option: " << vertex_identifier << "." << std::endl;
      std::terminate();
    }
  }
  // Closes the stream.
  output_file.close();
}


// ================================================================================================
// ================================================================================================
template<typename graph_t>
void notBGL::save_weighted_edgelist(std::string filename, graph_t &g, bool write_names = true)
{
  // Stream objects.
  std::ofstream edgelist_file;
  // Objects related to the weights.
  typename boost::property_map<graph_t, boost::edge_weight_t>::type weight_map = get(boost::edge_weight, g);
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
      for (std::tie(it, end) = boost::edges(g); it!=end; ++it)
      {
        edgelist_file << g[source(*it, g)].name << " "
                      << g[target(*it, g)].name << " "
                      << get(weight_map, *it)   << std::endl;
      }
    }
    else // vertices are identified using the numerical id.
    {
      for (std::tie(it, end) = edges(g); it!=end; ++it)
      {
        edgelist_file << g[source(*it, g)].id << " "
                      << g[target(*it, g)].id << " "
                      << get(weight_map, *it)   << std::endl;
      }
    }
  }
  // Closes the stream.
  edgelist_file.close();
}


// ================================================================================================
// ================================================================================================
template<typename graph_t, typename map_t>
auto notBGL::load_coordinates(std::string filename, graph_t &g, map_t &Name2Vertex, const unsigned int N = 3)
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
  std::vector<double> position(N);
  // Opens the stream and terminates if the operation did not succeed.
  positions_file.open(filename.c_str(), std::ios_base::in);
  if( !positions_file.is_open() )
  {
    std::cerr << "Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }
  else
  {
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
        vertex_t v1 = boost::add_vertex(g);
        Name2Vertex[name1_str] = v1;
        g[v1].name = name1_str;
        g[v1].id = boost::num_vertices(g);
        // std::cerr << "Vertex " << name1_str << " in the positions file is not"
        //           << " in the edgelist (has a degree equal to 0). The options"
        //           << " to solve this problem are"
        //           << std::endl << std::endl;
        // std::cerr << "  1) Remove the line corresponding to vertex "
        //           << name1_str << " from the file " << filename << "."
        //           << std::endl;
        // std::cerr << "  2) Add the line '" << name1_str << " " << name1_str
        //           << "' (a self-loop) at the end of the edgelist file. This"
        //           << " will works given that the function loading the graph"
        //           << " from the edgelist do not allow self-loops (e.g.,"
        //           << " notBGL::load_edgelist())."
        //           << std::endl;
        // std::cerr << "  3) Modify the code (where this message can be found)"
        //           << " to make it add the zero-degree vertex " << name1_str
        //           << " to the graph_t object."
        //           << std::endl << std::endl;
        // std::terminate();
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

// // ================================================================================================
// // ================================================================================================
// template<typename graph_t>
// auto notBGL::load_edgelist(std::string filename, graph_t &g)
// {
//   // Vertex descriptors.
//   typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
//   vertex_t v1, v2;
//   // Stream objects.
//   std::ifstream edgelist_file;
//   std::stringstream one_line;
//   // String objects.
//   std::string full_line, name1_str, name2_str;
//   // Integer objects.
//   unsigned int node_cnt(0);
//   // Map objects to assure the uniqueness of vertices.
//   std::map< std::string, vertex_t > Name2Vertex;
//   // Iterator objects.
//   typename std::map< std::string, vertex_t >::iterator name_it;
//   // Opens the stream and terminates if the operation did not succeed.
//   edgelist_file.open(filename.c_str(), std::ios_base::in);
//   if( !edgelist_file.is_open() )
//   {
//     std::cerr << "Could not open file: " << filename << "." << std::endl;
//     std::terminate();
//   }
//   else
//   {
//     // Reads the edgelist and registers the nodes and the edges.
//     while( !edgelist_file.eof() )
//     {
//       // Reads a line of the edgelist.
//       std::getline(edgelist_file,full_line); edgelist_file >> std::ws;
//       one_line.str(full_line); one_line >> std::ws;
//       one_line >> name1_str >> std::ws;
//       one_line >> name2_str >> std::ws;
//       one_line.clear();
//       // Is name1 new?
//       name_it = Name2Vertex.find(name1_str);
//       if( name_it == Name2Vertex.end() )
//       {
//         v1 = boost::add_vertex(g);
//         Name2Vertex[name1_str] = v1;
//         g[v1].name = name1_str;
//         g[v1].id = node_cnt;
//         ++node_cnt;
//       }
//       else
//       {
//         v1 = name_it->second;
//       }
//       // Is name2 new?
//       name_it = Name2Vertex.find(name2_str);
//       if( name_it == Name2Vertex.end() )
//       {
//         v2 = boost::add_vertex(g);
//         Name2Vertex[name2_str] = v2;
//         g[v2].name = name2_str;
//         g[v2].id = node_cnt;
//         ++node_cnt;
//       }
//       else
//       {
//         v2 = name_it->second;
//       }
//       // Creates the edge if it does not exist (forbids multiple edges and self-loops).
//       if(v1 != v2) // not a self-loop.
//       {
//         if(boost::edge(v1,v2,g).second == false) // not a multiple edge
//         {
//           boost::add_edge(v1, v2, g);
//         }
//       }
//     }
//   }
//   // Closes the stream.
//   edgelist_file.close();
//   // Returns the maps "Name2Vertex" to help add other properties to vertices.
//   return Name2Vertex;
// }