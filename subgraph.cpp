/**
 * @file subgraph.cpp
 * Test script for viewing a subgraph from our Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <iterator>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"

#include "Graph.hpp"

/** An iterator that skips over elements of another iterator based on whether
 * those elements satisfy a predicate.
 *
 * Given an iterator range [@a first, @a last) and a predicate @a pred,
 * this iterator models a filtered range such that all i with
 * @a first <= i < @a last and @a pred(*i) appear in order of the original range.
 */
template <typename Pred, typename It>
class filter_iterator
    : private equality_comparable<filter_iterator<Pred,It>> {
 public:
  // Get all of the iterator traits and make them our own
  typedef typename std::iterator_traits<It>::value_type        value_type;
  typedef typename std::iterator_traits<It>::pointer           pointer;
  typedef typename std::iterator_traits<It>::reference         reference;
  typedef typename std::iterator_traits<It>::difference_type   difference_type;
  typedef typename std::input_iterator_tag                     iterator_category;

  typedef filter_iterator<Pred,It> self_type;

  /** Constructor
  * @RI @a it_ != end_ && p_(*it_)
  */
  filter_iterator(const Pred& p, const It& first, const It& last)
      : p_(p), it_(first), end_(last) {
    // HW1 #4: YOUR CODE HERE
    fix();
  }

  // HW1 #4: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  value_type operator*() const
  {
    return (*it_);
  }
  
  self_type& operator++()
  {
    //std::cout << (*it_).index() << std::endl;
    ++it_;
    fix();
    return *this;
  }

  bool operator==(const self_type& x) const
  {
    return it_ == x.it_;
    //return ((it_ == x.it_) && (end_ == x.end_));
  }

 private:
  Pred p_;
  It it_;
  It end_;

  /*@pre Do not need to keep the RI 
  * @post fix current state so that meets RI
  */
  void fix()
  {
    while((it_ != end_) && (!p_(*it_)))
    {
      ++it_;
    }
  }
};

/** Helper function for constructing filter_iterators.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter>
filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
                                         const Pred& p) {
  return filter_iterator<Pred,Iter>(p, it, end);
}

/* Calculate the absolute value of a given number, original only works for integer*/
double myAbs(double x)
{
  return x < 0 ? -x : x;
}

// HW1 #4: YOUR CODE HERE

/* Predicate that makes a central cross of a given graph, which
*  slices a 3D model into eight parts.
*/
struct CrossPredicate {

  CrossPredicate(double r) : range_(r) {}

  template <typename NODE>
  bool operator()(const NODE& n) {
    bool xyPlane = (myAbs(n.position().z) < range_);
    bool xzPlane = (myAbs(n.position().y) < range_);
    bool yzPlane = (myAbs(n.position().x) < range_);
    return !(xyPlane || xzPlane || yzPlane);
  }
private:
  // specify inside which range nodes won't appear
  double range_;
};


/** Test predicate for HW1 #4 */
struct SlicePredicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return n.position().x < 0;
  }
};

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  typedef Graph<int> GraphType;
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // HW1 #4: YOUR CODE HERE
  // Use the filter_iterator to plot an induced subgraph.
  auto node_map = viewer.empty_node_map(graph);
  double r = 0.05;
  auto firstIt = make_filtered(graph.node_begin(), graph.node_end(), CrossPredicate(r));
  auto lastIt = make_filtered(graph.node_end(), graph.node_end(), CrossPredicate(r));
  viewer.add_nodes(firstIt, lastIt, node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();

  return 0;
}
