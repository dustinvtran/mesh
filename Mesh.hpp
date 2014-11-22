#pragma once
/** @file Mesh.hpp
 * @brief A Mesh is composed of nodes, edges, and triangles such that:
 *  -- All triangles have three nodes and three edges.
 *  -- All edges belong to at least one triangle and at most two triangles.
 */

#include <boost/iterator/transform_iterator.hpp>
#include <cmath>
#include "Graph.hpp"

/** @class Mesh
 * @brief A template for 3D triangular meshes.
 *
 * Users can add triangles and retrieve nodes, edges, and triangles.
 */
template <typename N, typename E, typename T>
class Mesh {
  // HW3: YOUR CODE HERE
  // Write all typedefs, public member functions, private member functions,
  //   inner classes, private member data, etc.
 public:
  /** Type of indexes and sizes. Return type of Mesh::num_nodes(). */
  typedef int size_type;
  
  typedef N node_value_type;
  typedef E graph_edge_value_type;
  typedef T triangle_value_type;
  
 private:
  struct internal_triangle;
  template <typename graph_edge_value_type>
  struct MeshEdgeValue {
    /* Indices of the left and right triangle of the edge. */
    size_type triangle_left;
    size_type triangle_right;

    graph_edge_value_type value;

    MeshEdgeValue() : triangle_left(-1), triangle_right(-1) {
      
    }
  };
  
  typedef Graph<N, MeshEdgeValue<graph_edge_value_type>> GraphType;
  
 public:
  typedef typename GraphType::Node Node;
  typedef Node node_type;
  
  typedef typename GraphType::Edge Edge;
  typedef Edge edge_type;
  
  typedef typename GraphType::incident_iterator incident_iterator;
  typedef typename GraphType::node_iterator node_iterator;
  typedef typename GraphType::edge_iterator edge_iterator;
  
  /* Iterator to iterator through all 3 adjacent triangles of a given triangle */
  class TriTriangleIterator;
  
  class Triangle : private totally_ordered<Triangle> {
   public:    
    /** Access the ith node of the triangle.
      * @param[in] i 0, 1, or 2
      *
      * Complexity: O(1).
      */
    Node node(size_type i) const {
      assert(0 <= i && i <= 2);
      return mesh_->graph_.node(fetch().nodes[i]);
    }
    
    /** Access the ith edge of the triangle.
      * @param[in] i 0, 1, or 2
      *
      * Complexity: O(1).
      */
    Edge edge(size_type i) const {
      assert(0 <= i && i <= 2);
      if (i == 0) {
        return mesh_->graph_.edge(node(0), node(1));
      } else if (i == 1) {
        return mesh_->graph_.edge(node(0), node(2));
      } else {
        return mesh_->graph_.edge(node(1), node(2));
      }
    }

    /** Calculate the area of the triangle.
      *
      * Complexity: O(1).
      */
    double area() const {
      // using heron's formula
      double s = (edge(0).length() + edge(1).length() + edge(2).length())/2;
      return std::sqrt(s *
                      (s - edge(0).length()) *
                      (s - edge(1).length()) *
                      (s - edge(2).length())
                      );
    }

    /** Access the index of the triangle.
      *
      * Complexity: O(1).
      */
    size_type index() const {
      return index_;
    }

    /** Access the value(s) stored in the triangle.
      *
      * Complexity: O(1).
      */
    triangle_value_type& value() {
      return fetch().value;
    }

    /** Access the value(s) stored in the triangle.
      *
      * Complexity: O(1).
      */
    const triangle_value_type& value() const {
      return fetch().value;
    }
    
    bool operator==(const Triangle& x) const {
      return std::tie(mesh_, node(0), node(1), node(2)) == std::tie(x.mesh_, x.node(0), x.node(1), x.node(2));
    }
    
    bool operator<(const Triangle& x) const {
      return std::tie(mesh_, node(0), node(1), node(2)) < std::tie(x.mesh_, x.node(0), x.node(1), x.node(2));
    }


    /** Return the iterator which corresponds to the first adjacent triangle. */
    TriTriangleIterator begin() const {
      return {mesh_, index_, 0};
    }

    /** Return the iterator which corresponds to the past-the-last adjacent triangle. */
    TriTriangleIterator end() const {
      return {mesh_, index_, 3};
    }

    /** Construct invalid triangle. */
    Triangle() : mesh_(nullptr), index_(-1) {}
    
    bool valid() {
      return index_ != -1;
    }

  private:
    friend class Mesh;
    /** Private valid constructor */
    Triangle(const Mesh* m, size_type i) : mesh_(const_cast<Mesh*>(m)), index_(i) {
    }
    /** Pointer to the Mesh object this triangle belongs */
    Mesh* mesh_;
    /** Index of this triangle, we do not distinguish UID from index in Triangle */
    size_type index_;

    internal_triangle& fetch() const {
      return mesh_->triangles_[index_];
    }
  };
  
  /** Return the number of nodes in the mesh. */
  size_type num_nodes() const {
    return graph_.num_nodes();
  }

  /** Return the number of edges in the mesh. */
  size_type num_edges() const {
    return graph_.num_edges();
  }

  /** Return the number of triangles in the mesh. */
  size_type num_triangles() const {
    return triangles_.size();
  }
  
  /** Add a triangle.
   * @param[in] n1, n2, n3 Valid nodes belonging to @a graph_.
   * @pre n1.index() != n2.index() != n3.index()
   * @pre For all 1 <= i < j <= 3, @a graph_.edge(ni, nj) has at most 1 triangle; or if graph_.edge(ni, nj) has 2 triangles, then has_triangle(ni, nj, nk)
   * @post If !has_triangle(n1, n2, n3),
   *   new @a num_triangles() == old @a num_triangles() + 1
   *   For all 1 <= i < j <= 3, @a graph_.edge(ni, nj) has 1 more triangle
   * Else
   *   No members are modified.
   * @post For all 1 <= i < j < k <= 3, if nk is on the left(right) side of the hyperplane formed by the vector from ni.position() to nj.position(), where ni < nj,
   *   then graph_.edge(ni, nj).left(right) == triangle(ni, nj, nk)
   * @return Triangle with vertices n1, n2, n3.
   *
   * Complexity: O(1).
   */
  Triangle add_triangle(const Node& n1, const Node& n2, const Node& n3, const triangle_value_type& val = triangle_value_type()) {
    /*
    1. Create a new internal_triangle object with n1, n2, n3,
        which will be pushed into @a triangles_;
    2. Add edges between each combination of the 3 nodes into @a graph_;
    3. For each edge_ij (i < j, total of 3 edges, e.g: A->B, A->C, B->C)
      a. calculate if the third node_k is on the left or right side towards
          this edge;
      b. set edge_ij.value().left(right)_tri to the index of the new triangle
      c. set edge_ij.value().left(right)_node to the index of the third node
    */
    assert(n1 != n2);
    assert(n1 != n3);
    assert(n2 != n3);
    if (has_triangle(n1, n2, n3)) {
      return triangle(n1, n2, n3);
    }
    // Pre-process the nodes, sort them so that the relationship if temp[0] < temp[1] < temp[2]
    auto new_triangle = internal_triangle(val);
    new_triangle.nodes.push_back(n1.index());
    new_triangle.nodes.push_back(n2.index());
    new_triangle.nodes.push_back(n3.index());
    std::sort(new_triangle.nodes.begin(), new_triangle.nodes.end());
    new_triangle.index = triangles_.size();
    triangles_.push_back(new_triangle);
    // Add edge (if it doesn't already exist), and assign the triangle toward whichever
    // side of the edge it faces.
    for (unsigned i=0; i < 3; ++i) {
        auto edge = graph_.add_edge(node(new_triangle.nodes[i % 3]), node(new_triangle.nodes[(i+1) % 3]));
        if (on_right(edge, node(new_triangle.nodes[(i+2) % 3]).position())) {
          assert(edge.value().triangle_right == -1);
          edge.value().triangle_right = new_triangle.index;
        } else {
          assert(edge.value().triangle_left == -1);
          edge.value().triangle_left = new_triangle.index;
       }
    }

    return triangle(new_triangle.index);
  }
  
  /** Verify whether there is a triangle with vertices n1, n2, n3.
    *
    * Complexity: O(1).
    */
  bool has_triangle(const Node& n1, const Node& n2, const Node& n3) const {
      if (!(graph_.has_edge(n1, n2)) || !(graph_.has_edge(n1, n3)) || !(graph_.has_edge(n2, n3))) {
        return false;
      }

      Edge edge = graph_.edge(n1, n2);
      auto tri_left = triangle_left(edge);
      auto tri_right = triangle_right(edge);
      if (tri_left.valid()) {
        return ((n3 == tri_left.node(0) || n3 == tri_left.node(1) || n3 == tri_left.node(2)));
      } else if (tri_right.valid()) {
        return ((n3 == tri_right.node(0) || n3 == tri_right.node(1) || n3 == tri_right.node(2)));
      }
      return false;
  }

  /** Access the triangle with vertices n1, n2, n3.
    * @pre has_triangle(n1, n2, n3)
    *
    * Complexity: O(1).
    */
  Triangle triangle(const Node& n1, const Node& n2, const Node& n3) const {
    /*
    Pick up an arbitrary edge of the three nodes, check if this edge's
      left(right)_node is equal to the index of the third node
    */
    if (!has_triangle(n1, n2, n3)) {
        return Triangle();
    }    
    Edge edge = graph_.edge(n1, n2);
    auto tri_left = triangle_left(edge);
    // check if the third node @a n3 is on the left triangle of the edge
    if (tri_left.valid() && 
        (n3 == tri_left.node(0) || n3 == tri_left.node(1) || n3 == tri_left.node(2))) {
        return tri_left;
    } else {
        return triangle_right(edge);
    }
  }
  
  Triangle triangle(size_type i) const {
    return {this, i};
  }
  
  /** Add a node.
   * Satisfies the conditions for @a graph_'s add_node().
   *
   * Complexity: O(1).
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    return graph_.add_node(position, val);
  }

  Node node(size_type i) {
    return graph_.node(i);
  }
  
  /** Return the left triangle of an edge.
    *
    * Complexity: O(1).
    */
  Triangle triangle_left(const Edge& edge) const {
    if (edge.value().triangle_left == -1) {
      return Triangle();
    }
    
    return triangle(edge.value().triangle_left);
  }

  /** Return the right triangle of an edge.
    *
    * Complexity: O(1).
    */
  Triangle triangle_right(const Edge& edge) const {
    if (edge.value().triangle_right == -1) {
      return Triangle();
    }
    
    return triangle(edge.value().triangle_right);
  }
  
  /** Return the outward normal of an edge given the triangle.
    * @param[in] triangle A valid triangle containing @a edge in the mesh.
    * @param[in] edge A valid edge.
    * @return Point which specifies the outward normal vector.
    *
    * Complexity: O(1).
    */
  Point out_norm(const Triangle& triangle, const Edge& edge) const {
    double dx = edge.node2().position().x - edge.node1().position().x;
    double dy = edge.node2().position().y - edge.node1().position().y;
    // TODO: normalize, and also maybe flip them
    double scale = sqrt(dx * dx + dy * dy);
    dx = dx / scale * edge.length();
    dy = dy / scale * edge.length();
    
    if (on_right(edge, Point(edge.node1().position().x-dy, edge.node1().position().y+dx, 0))
     && (triangle.index() == edge.value().triangle_left)) {
      return Point(-dy, dx, 0);
    } else {
      return Point(dy, -dx, 0);
    }
  }
#if 0
  /** Functor to use Boost::transform_iterator */
  struct InternalToTriangleTransform {
    Triangle operator()(internal_triangle internal_tri) {
      return mesh_->triangle(internal_tri.index);
    }
  private:
    friend class Mesh;
    InternalToTriangleTransform(Mesh* m) : mesh_(m) {
    }

    Mesh* mesh_;
  };

  /** Define the TriangleIterator type */
  using TriangleIterator = boost::transform_iterator<InternalToTriangleTransform, typename std::vector<internal_triangle>::const_iterator>;
  
  TriangleIterator triangle_begin() const {
      return {triangles_.begin(), InternalToTriangleTransform(this)};
  }
  
  TriangleIterator triangle_end() const {
      return {triangles_.end(), InternalToTriangleTransform(this)};
  }
#else
  class TriangleIterator : private totally_ordered<TriangleIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct invalid iterator. */
    TriangleIterator() : mesh_(nullptr), index_(-1) { 
    
    }

    /** Return the triangle at the iterator. */
    Triangle operator*() {
      return mesh_->triangle(index_);
    }

    /** Increment iterator. */
    TriangleIterator operator++() {
      if (index_ < mesh_->num_triangles())
        ++index_;
      return(*this);
    }

    /** Compare if two iterators are the same. */
    bool operator==(const TriangleIterator& x) const {
      return std::tie(mesh_, index_) == std::tie(x.mesh_, x.index_);
    }

  private:
    friend class Mesh;
    
    TriangleIterator(const Mesh* m, size_type idx) : 
      mesh_(const_cast<Mesh*>(m)), index_(idx) {
    }

    Mesh* mesh_;
    size_type index_;
  };

  /** Return the iterator which corresponds to the first triangle. */
  TriangleIterator triangle_begin() {
    return {this, 0};
  }

  /** Return the iterator which corresponds to the past-the-last triangle. */
  TriangleIterator triangle_end() {
    return {this, num_triangles()};
  }
#endif
  node_iterator node_begin() {
    return graph_.node_begin();
  }
  
  node_iterator node_end() {
    return graph_.node_end();
  }
  
  edge_iterator edge_begin() {
    return graph_.edge_begin();
  }
  
  edge_iterator edge_end() {
    return graph_.edge_end();
  }

  class NodeTriangleIterator : private totally_ordered<NodeTriangleIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct invalid iterator. */
    NodeTriangleIterator() :
      mesh_(nullptr), it_(incident_iterator()), it_end_(incident_iterator()) {
    }

    /** Return the adjacent triangle at the iterator. */
    Triangle operator*() {
      auto edge = (*it_);
        size_type tri_index = edge.node1() < edge.node2() ?
                 mesh_->triangle_left(edge).index() : mesh_->triangle_right(edge).index();
      return mesh_->triangle(tri_index);
    }

    /** Increment iterator. */
    NodeTriangleIterator operator++() {
      if (it_ != it_end_) {
        ++it_;
        fix();
      }
      return (*this);
    }

    /** Compare if two iterators are the same. */
    bool operator==(const NodeTriangleIterator& x) const {
      return std::tie(mesh_, it_, it_end_) == std::tie(x.mesh_, x.it_, it_end_);
    }

   private:
    friend class Mesh;
    NodeTriangleIterator(const Mesh* m, incident_iterator it, incident_iterator it_end) :
      mesh_(const_cast<Mesh*>(m)), it_(it), it_end_(it_end) {
      fix();
    }

    void fix() {
      while (it_ != it_end_) {
        auto edge = (*it_);
        size_type tri_index = edge.node1() < edge.node2() ?
                 mesh_->triangle_left(edge).index() : mesh_->triangle_right(edge).index();
        if (tri_index == -1) {
          ++it_;
        } 
        else {
          break;
        }
      }
    }

    Mesh* mesh_;
    incident_iterator it_;
    incident_iterator it_end_;
  };

  /** Return the iterator which corresponds to the first adjacent triangle. */
  NodeTriangleIterator adjacent_triangle_begin(const Node& node) const {
    return {this, node.edge_begin(), node.edge_end()};
  }

  /** Return the iterator which corresponds to the past-the-last adjacent triangle. */
  NodeTriangleIterator adjacent_triangle_end(const Node& node) const {
    return {this, node.edge_end(), node.edge_end()};
  }

  class TriTriangleIterator : private totally_ordered<TriTriangleIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;


    /** Construct invalid iterator. */
    TriTriangleIterator() :
      mesh_(nullptr), triangle_index_(-1), edge_index_(-1) {
    }

    /** Return the adjacent triangle at the iterator. */
    Triangle operator*() {
      Edge edge = mesh_->triangle(triangle_index_).edge(edge_index_);
      return mesh_->triangle_left(edge).index() == triangle_index_ ?
        triangle_right(edge) : triangle_left(edge);
    }

    /** Increment iterator. */
    TriTriangleIterator operator++() {
      if (edge_index_ < 3)
        ++edge_index_;
      return(*this);
    }

    /** Compare if two iterators are the same. */
    bool operator==(const TriTriangleIterator& x) const {
      return std::tie(mesh_, triangle_index_, edge_index_) ==
             std::tie(x.mesh_, x.triangle_index_, x.edge_index_);
    }

  private:
    friend class Mesh;
    TriTriangleIterator(const Mesh* m, size_type tri_idx, size_type e) :
      mesh_(const_cast<Mesh*>(m)), triangle_index_(tri_idx), edge_index_(e) {
    }

    Mesh* mesh_;
    size_type triangle_index_;
    size_type edge_index_;
  };  

 private:
  /** Helper function to check if @a pt is on the right side of @a edge .
    * The direction of the @a edge is defined as from min(edge.node1(), edge.node2()) to max(edge.node1(), edge.node2())
    *
    * @pre @a pt cannot be on the line of @a edge
    */
  bool on_right(const Edge& edge, const Point& pt) const {
    auto A_pos = (std::min(edge.node1(), edge.node2())).position();
    auto B_pos = (std::max(edge.node1(), edge.node2())).position();
    double calculated_side = (B_pos.x - A_pos.x) * (pt.y - A_pos.y) - 
                          (B_pos.y - A_pos.y) * (pt.x - A_pos.x);

    assert(calculated_side != 0);
    return 0 < calculated_side;
  }

  GraphType graph_;
  /* Internal data structure for triangle*/
  struct internal_triangle {  
    internal_triangle(const triangle_value_type& val) : value(val) {
    }
    
    /* Vector to stores three indices of the nodes on the triangle */
    std::vector<size_type> nodes;
    /** Index of the triangle */
    size_type index;
    triangle_value_type value;
  };
  // can be used for triangle traversal
  std::vector<internal_triangle> triangles_;
};