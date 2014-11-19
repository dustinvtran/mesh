/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SDLViewer to visualize the solution.
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"

#include "Graph.hpp"
#include "Point.hpp"
#include "BoundingBox.hpp"

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

typedef Graph<double, char> GraphType;
typedef double value_type;

/* Boundary type 1 object
* Time complexity: O(1)
*/
struct Bound1 {
  unsigned operator()(const Point& x) {
    return (norm_inf(x) == 1) ? 1 : 0;
  }
};

/* Boundary type 1 object
* Time complexity: O(1)
*/
struct Bound2
{
  unsigned operator()(const Point& x) {
    Point pt1 = Point(0.6, 0.6, 0);
    Point pt2 = Point(-0.6, 0.6, 0);
    bool on_bound = ((norm_inf(x - pt1) < 0.2) || (norm_inf(x + pt1) < 0.2) ||
          (norm_inf(x - pt2) < 0.2) || (norm_inf(x + pt2) < 0.2));

    return on_bound ? 2 : 0;
  }
};

/* Boundary type 1 object
* Time complexity: O(1)
*/
struct Bound3
{
  unsigned operator()(const Point& x) {
    BoundingBox cond3Bb = BoundingBox(Point(-0.6,-0.2,-1), Point(0.6,0.2,1));
    return cond3Bb.contains(x) ? 3 : 0;
  }
};

/* Struct to combine different boundary condition together, each Boundary condition
*   should return an unsigned index, which represents which condition index if the 
*   node is meeting or not.
* @param[in] x A Point object
* @return 0: @a x does not meet any of the condition; otherwise return a none zero unsigned.
*
* @pre All boundary condition functor should accept a Point object and return an 
*      unsigned value, where > 0 means the Point meets this boundary condition, == 0 otherwise
* @pre The intersection between any two boundary condition should be an empty set.
* 
* @post Time Complexity(TC): O(max(TC(b1), TC(b2)))
*/
template <typename B1, typename B2>
struct MetaBound
{
  MetaBound(B1 b1, B2 b2) : b1_(b1), b2_(b2) {

  }

  unsigned operator()(const Point& x) {
    return std::max(b1_(x), b2_(x));
  }
private:
  B1 b1_;
  B2 b2_;
};

/* Syntax sugar function for making a MetaBound object from two object satisfies Boundary concept*/
template <typename B1, typename B2>
MetaBound<B1, B2> make_combined_boundary(B1 b1, B2 b2) {
  return {b1, b2};
}
/* Function to compute f(.) value
* @param[in] x A Point object
*/
value_type f(const Point& x) {
  value_type fx = 5 * cos(norm_1(x));
  return fx;
}

/* Function to compute g(.) value
* @param[in] x A Point object
* @param[in] bd A bound object to check which boundary condition @a x satisfies
*
* @pre @a x satisfies and only satisfies one of the boundary conditions
*/
template <typename Bound>
value_type g(const Point& x, Bound& bd) {
  value_type bound_scores[4] = {0, 0, -0.2, 1};
  unsigned bound_type = bd(x);

  assert(bound_type > 0 && bound_type < 4);
  return bound_scores[bound_type];
}

/* Self defined Matrix class to compute A*x for HW3 problem 
* It is a Has-a relationship with GraphType
*/
class GraphSymmetricMatrix {
public:
  GraphSymmetricMatrix(GraphType& g) : graph_(g) {
  }

  /** Helper function to perform multiplication . Allows for delayed
  * evaluation of results .
  * Assign :: apply (a , b ) resolves to an assignment operation such as
  * a += b , a -= b , or a = b .
  * @pre @a size ( v ) == size ( w ) 
  * 
  * I don't know the real mechanism. But if I am designing this, I will reload operator=, +=, -=
  * defined in mtl::vec::mat_cvec_multiplier<Mat, Vec> so that when such operators are called,
  * mat_cvec_multiplier will call A.mult(v, w, proper_Assign) inside.
  */
  template <typename VectorIn, typename VectorOut, typename Assign>
  void mult(const VectorIn& v, VectorOut& w, Assign) const {
    assert(std::size_t(size(v)) == mat_size());
    assert(size(v) == size(w));
    typedef std::size_t size_type;
    
    auto bounds = make_combined_boundary(make_combined_boundary(Bound1(), Bound2()), Bound3());
    
    for (auto node_i : graph_.nodes()) {
      value_type sum = 0;
      /* if node_i is on boundary, as long as j!=i, A(i,j) = 0,
      * even if j==i, A(i,j)=1, so sum definately is equal to 1 * v[i]
      * Time complexity: O(1)
      */
      if (bounds(node_i.position())) {
        sum = v[node_i.index()];
      }
      /* else: node_i is not on boundary
      * For this case, we only need to check every incident node to node_i because otherwise, A_ij = 0
      * Time complexity: O(degree)
      */
      else {
        // since node_i is not on boundary, A_ii = -degree(i)
        sum -= node_i.degree() * v[node_i.index()];
        for (auto icdtEdge : node_i.incidentEdges()) {
          auto node_j = icdtEdge.node2();
          // We should only add v[j] only when node_j is not on boundary since otherwise, A_ij = 0.
          // If both i and j are not on the boundary && edge(i,j) exists, A_ij = 1
          if (!bounds(node_j.position())) {
            sum += v[node_j.index()];
          }
        }
      }
      Assign::apply(w[node_i.index()], sum);
    }
  }

  /* * Matrix - vector multiplication forwards to MTL ’s lazy mat_cvec_multiplier */
  template <typename VectorIn>
  mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, VectorIn> operator*(const VectorIn& v) const {
    return {*this, v};
  }

  std::size_t mat_size() const {
    return graph_.size();
  }

private:
  GraphType& graph_;
};

inline std::size_t size(const GraphSymmetricMatrix& A) {
  return A.mat_size() * A.mat_size();
}

inline std::size_t num_rows(const GraphSymmetricMatrix& A) {
  return A.mat_size();
}

inline std::size_t num_cols(const GraphSymmetricMatrix& A) {
  return A.mat_size();
}

/* * Traits that MTL uses to determine properties of GraphSymmetricMatrix. */
namespace mtl {
  /** GraphSymmetricMatrix implements the Collection concept
  *   with value_type and size_type */
  template <>
  struct Collection<GraphSymmetricMatrix>
  {
    typedef double value_type;
    typedef int size_type;
  };

  namespace ashape {
  /** Define IdentityMatrix to be a non - scalar type . */
    template<>
    struct ashape_aux<GraphSymmetricMatrix>
    {       
      typedef nonscal type;
    };
  }
}

// Color functor that returns a CS207::Color according to a node's value
struct HeatmapColor {
  /* Calculate the head color
  * @pre @a min_ <= @a n.value() <= @a max_
  */
  template <typename NODE>
  CS207::Color operator()(const NODE& n) {
    if (max_ == min_ || n.value() == min_)
      return CS207::Color::make_heat(0);
    value_type val = (n.value() - min_) / (max_ - min_);
    // in case that some nodes' value is still -1, normalized it.
    return CS207::Color::make_heat(val);
  }

  HeatmapColor(value_type min, value_type max) : min_(min), max_(max) {
    assert(min_ <= max_);
  }
private:
  value_type min_;
  value_type max_;
};

/* Functor to return the position of a node where its original z coordinate
* is replaced by node.value() */
struct ValuePosition {
  template <typename NODE>
  Point operator()(const NODE& node) {
    Point pos = node.position();
    pos.z = node.value();
    return pos;
  }
};

/* Self defined iteration class to visualize Graph during conjugate descending*/
template <typename Viewer, typename NodeMap, typename Vector, class Real, class OStream = std::ostream>
class visual_iteration : public itl::cyclic_iteration<Real, OStream> {
  typedef itl::cyclic_iteration<Real, OStream> super;
  typedef visual_iteration self;

private:
  void visualize_result(const Vector& v) const {
    assert(size(v) == graph_.size());

    viewer_.clear();
    node_map_.clear();

    for (auto node : graph_.nodes()) {
      node.value() = v[node.index()];
    }

    value_type minValue = *std::min_element(v.begin(), v.end());
    value_type maxValue = *std::max_element(v.begin(), v.end());

    node_map_ = viewer_.empty_node_map(graph_);
    viewer_.add_nodes(graph_.node_begin(), graph_.node_end(), HeatmapColor(minValue, maxValue), ValuePosition(), node_map_);
    viewer_.add_edges(graph_.edge_begin(), graph_.edge_end(), node_map_);
  }
  /* Retained graph reference */
  GraphType& graph_;
  /* Retained viewer reference */
  Viewer& viewer_;
  /* Retained node map reference */
  NodeMap& node_map_;
  /* Retained vector reference, which is being solved */
  Vector& sol_;

public:
  /* Constructor */
  visual_iteration(GraphType& graph, Viewer& viewer, NodeMap& node_map, 
                  const Vector& r0, Vector& sol, int max_iter_, Real tol_, 
                  Real atol_ = Real(0), int cycle_ = 100, OStream& out = std::cout)
  : super(r0, max_iter_, tol_, atol_, cycle_, out),
  graph_(graph), viewer_(viewer), node_map_(node_map), sol_(sol) {

  }
  // Personally, I don't think we need to overwrite this function
  bool finished() {
    return super::finished();
  }

  template <typename T>
  bool finished(const T& v) {
    bool ret = super::finished(v);
    visualize_result(sol_);
    return ret;
  }
};

/* Function to make a visualize_iteration object */
template <typename Viewer, typename NodeMap, typename Vector, class Real, class OStream = std::ostream>
visual_iteration<Viewer, NodeMap, Vector, Real, OStream> make_visualize_iterator (
  const Vector& r0, Vector& sol, int max_iter_, Real tol_, Real atol_, int cycle_,
  GraphType& graph, Viewer& viewer, NodeMap& node_map,  OStream& out = std::cout) {
  return {graph, viewer, node_map, r0, sol, max_iter_, tol_, atol_, cycle_, out};
}

/** Remove all the nodes in graph @a g whose posiiton is contained within
 * BoundingBox @a bb
 * @post For all i, 0 <= i < @a g.num_nodes(), !bb.contains(g.node(i).position())
 */
void remove_box(GraphType& g, const BoundingBox& bb) {
  auto it = g.node_begin();
  while(it != g.node_end()) {
    if (bb.contains((*it).position()))
      it = g.remove_node(it);
    else
      ++it;
  }
}

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE EDGES_FILE\n";
    exit(1);
  }

  // Define an empty Graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<typename GraphType::node_type> node_vec;
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    node_vec.push_back(graph.add_node(2*p - Point(1,1,0)));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t)) {
    graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
    graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
    graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
    graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
  }

  // Get the edge length, should be the same for each edge
  double h = graph.edge(0).length();

  // Make holes in our Graph
  remove_box(graph, BoundingBox(Point(-0.8+h,-0.8+h,-1), Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, BoundingBox(Point( 0.4+h,-0.8+h,-1), Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, BoundingBox(Point(-0.8+h, 0.4+h,-1), Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, BoundingBox(Point( 0.4+h, 0.4+h,-1), Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, BoundingBox(Point(-0.6+h,-0.2+h,-1), Point( 0.6-h, 0.2-h,1)));

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;
  
  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  auto node_map = viewer.empty_node_map(graph);
  // Define b using the graph, f, and g.
  // Construct the GraphSymmetricMatrix A using the graph
  // Solve Au = b using MTL.

  GraphSymmetricMatrix A = GraphSymmetricMatrix(graph);
  itl::pc::identity<GraphSymmetricMatrix> P(A);

  mtl::dense_vector<value_type> x(A.mat_size()), b(A.mat_size());

  auto bounds = make_combined_boundary(make_combined_boundary(Bound1(), Bound2()), Bound3());

  for (unsigned long i = 0; i < size(b); ++i) {
    if (bounds(graph.node(i).position())) {
      b[i] = g(graph.node(i).position(), bounds);
    }
    else {
      b[i] = h * h * f(graph.node(i).position());
      for (auto icdtEdge : graph.node(i).incidentEdges()) {
        auto icdtNode = icdtEdge.node2();
        if (bounds(graph.node(icdtNode.index()).position())) {
          b[i] -= g(icdtNode.position(), bounds);
        }
      }
    }
  }

  //b = A * x;
  //itl::cyclic_iteration<value_type> iter(b, 1e3, 1.e-10, 0.0, 50);
  auto visual_iter = make_visualize_iterator(b, x, 4e2, 1.e-10, 0.0, 50, graph, viewer, node_map);
  cg(A, x, b, P, visual_iter);

  value_type minValue = *std::min_element(x.begin(), x.end());
  value_type maxValue = *std::max_element(x.begin(), x.end());

  viewer.add_nodes(graph.node_begin(), graph.node_end(), HeatmapColor(minValue, maxValue), ValuePosition(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  return 0;
}
