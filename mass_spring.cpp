/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point velocity;  //< Node velocity
  double mass;     //< Node mass
  // Damping coefficient
  double c;

  NodeData() : velocity(Point(0,0,0)), mass(0), c(1.0) {

  }
};

struct EdgeData
{
  // Spring coefficient
  double K;
  // Rest length of edge
  double L;

  EdgeData() : K(100), L(0) {

  }
};

// HW2 #1 YOUR CODE HERE
// Define your Graph type
typedef Graph<NodeData, EdgeData> GraphType;
typedef typename GraphType::node_type Node;
typedef typename GraphType::edge_type Edge;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on Node
 *           at time @a t.
 */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the {n+1} node positions
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().velocity * dt;
  }

  // Compute the {n+1} node velocities
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().velocity += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force being applied to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  Point operator()(Node n, double) {
    // HW2 #1: YOUR CODE HERE
    if (n.position().y == 0 && n.position().z == 0 && 
      (n.position().x == 1 || n.position().x == 0)) 
      return Point(0,0,0);

    //double K = 100.0;

    Point force = Point(0, 0, -grav) * n.value().mass;
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      Point nd2Pos = (*it).node2().position();
      //std:: cout << "K: " << (*it).value().K << std::endl;
      //std::cout << "rest L: " << (*it).value().L << std::endl;
      Point ijForce = -(*it).value().K * (n.position() - nd2Pos) * ((*it).length() - (*it).value().L) / (*it).length();
      force += ijForce;
    }
    
    return force;
  }
};

/** Force function object for Gravity. */
struct GravityForce
{
  Point operator()(Node n, double t) {
    (void)t;
    return Point(0, 0, -grav) * n.value().mass;
  }
};

/** Force function object for Mass Spring. */
struct MassSpringForce
{
  Point operator()(Node n, double t) {
    (void)t;
    Point force = Point(0, 0, 0);
    // for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
    for (auto edge : n.incidentEdges()) {
      Point nd2Pos = edge.node2().position();
      Point ijForce = -edge.value().K * (n.position() - nd2Pos) * (edge.length() - edge.value().L) / edge.length();
      force += ijForce;
    }
    return force;
  }
};

/** Force function object for damping force. */
struct DampingForce
{
  Point operator()(Node n, double t) {
    (void)t;
    return -n.value().velocity * n.value().c;
  }
};

/** Force function object for a special wind force. */
struct WindForce {
  Point operator()(Node n, double t) {
    Point force = Point(0, 5*t*n.position().x, 0);
    return force * n.value().mass;
  }
};

/** Force function object for combining two Force functor. */
template <typename F1, typename F2>
struct MetaForce
{
  MetaForce(F1 f1, F2 f2) : f1_(f1), f2_(f2) {
  }

  Point operator()(Node n, double t) {
    return f1_(n, t) + f2_(n, t);
  }
private:
  F1 f1_;
  F2 f2_;
};

/* Combine the effects of two forces. A Force is a struct that implements operator(),
* which takes in a Node and double (stands for time), and returns a Point
* 
* @param[in] f1, f2 Two force object
* @return A Metaforce object whose output force will be the summation of @a f1 and @f2
* @type F1 and F2 should implment Force functor's concept.
*/
template <typename F1, typename F2>
MetaForce<F1, F2> make_combined_force(F1 f1, F2 f2) {
  return MetaForce<F1, F2>(f1, f2);
}

/* Combine the effects of three forces.
* 
* @param[in] f1, f2, f3 Two force object
* @return A Metaforce object whose output force will be the summation of three forces.
* @type F1 and F2 should implment Force functor's concept.
*/
template <typename F1, typename F2, typename F3>
MetaForce<MetaForce<F1, F2>, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  return make_combined_force(make_combined_force(f1, f2), f3);
}

/** Constraint on Point(0,0,0) and Point(1,0,0) by fixing them */
template <typename G>
struct  PointConstraint
{
  void operator()(G& g, double t) {
    (void)t;
    for (auto node : g.nodes()) {
      if (norm(node.position() - Point(0,0,0)) < 0.01) {
        node.position() = Point(0,0,0);
        node.value().velocity = Point(0,0,0);
      }
      else if (norm(node.position() - Point(1,0,0)) < 0.01) {
        node.position() = Point(1,0,0);
        node.value().velocity = Point(0,0,0);
      }
    }
  }
};

/** Constraint on Plane z=-0.75, no point shall pass this plane */
template <typename G>
struct PlaneConstraint
{
  void operator()(G& g, double t) {
    (void)t;
    for (auto node : g.nodes()) {
      if (dot(node.position(), Point(0,0,1)) < z_) {
        Point new_pos = node.position();
        new_pos.z = z_;
        node.position() = new_pos;

        Point new_vel = node.value().velocity;
        new_vel.z = 0;
        node.value().velocity = new_vel;
      }
    }
  }
private:
  static double z_;
};

template <typename G>
double PlaneConstraint<G>::z_ = -0.75;

/** Constraint of a sphere, no point shall pass the sphere */
template <typename G>
struct SphereConstraint
{
  void operator()(G& g, double t) {
    (void)t;
    for (auto node : g.nodes()) {
      if (norm(node.position() - center_) < radius_) {
        Point r_i = (node.position() - center_) / norm(node.position() - center_);
        node.position() = center_ + r_i * radius_;

        Point old_vel = node.value().velocity;
        node.value().velocity = old_vel - dot(old_vel, r_i) * r_i;
      }
    }
  }
private:
  static double radius_;
  static Point center_;
};

template <typename G>
double SphereConstraint<G>::radius_ = 0.15;
template <typename G>
Point SphereConstraint<G>::center_ = Point(0.5,0.5,-0.5);

/** Constraint of a sphere, any point passes it shall be removed (burned) */
template <typename G>
struct SphereBurnConstraint
{
  void operator()(G& g, double t) {
    (void)t;
    // Cannot use Range because erase() requires an iterator!
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if (norm((*it).position() - center_) < radius_ * (t+0.5)) {
        it = g.remove_node(it);
      }
    }
  }
private:
  static double radius_;
  static Point center_;
};

template <typename G>
double SphereBurnConstraint<G>::radius_ = 0.15;
template <typename G>
Point SphereBurnConstraint<G>::center_ = Point(0.5,0.5,-0.5);

/** Constraint functor for combining two Constraint functor. */
template <typename C1, typename C2, typename G = GraphType>
struct MetaConstraint
{
  MetaConstraint(C1 c1, C2 c2) : c1_(c1), c2_(c2) {
  }

  void operator()(G& g, double t) {
    c1_(g, t);
    c2_(g, t);
  }
private:
  C1 c1_;
  C2 c2_;
};

/* Combine the effects of two constraints. A Constraint is a struct that implements operator(),
* which takes in a Graph&, modifies some internal data that violates constraints 
* while return nothing.
* 
* @param[in] c1, c2 Two constraint object
* @return A MetaConstraint object whose constraint effect will be the combination of @a c1 and @c2
* @type C1 and C2 should implment Constraint functor's concept.
* @pre The constraints @a c1 and @a c2 which will put onto Graph do not conflict with each other.
*/
template <typename C1, typename C2, typename G = GraphType>
MetaConstraint<C1, C2, G> make_combined_constraint(C1 c1, C2 c2) {
  return MetaConstraint<C1, C2, G>(c1, c2);
}

/* Combine the effects of three constraints.
* 
* @param[in] c1, c2, c3 Three constraint object.
* @return A MetaConstraint object whose constraint effect will be the combination of three constraints.
* @type C1, C2 and C3 should implment Constraint functor's concept.
* @pre The three constraints which will put onto Graph do not conflict with each other.
*/
template <typename C1, typename C2, typename C3, typename G = GraphType>
MetaConstraint<MetaConstraint<C1, C2, G>, C3, G> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
  return make_combined_constraint(make_combined_constraint(c1, c2), c3);
}

int main(int argc, char** argv) {
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<Node> nodes;
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t)) {
    for (unsigned i = 1; i < t.size(); ++i) {
      graph.add_edge(nodes[t[0]], nodes[t[1]]);
      graph.add_edge(nodes[t[0]], nodes[t[2]]);
#if 1
      // Diagonal edges: include as of HW2 #2
      graph.add_edge(nodes[t[0]], nodes[t[3]]);
      graph.add_edge(nodes[t[1]], nodes[t[2]]);
#endif
      graph.add_edge(nodes[t[1]], nodes[t[3]]);
      graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }
  }

  // Set initial conditions for nodes and edges.
  for (auto edge : graph.edges()) {
    edge.value().K = 100;
    edge.value().L = edge.length();
  }

  for (auto node : graph.nodes()) {
    node.value().mass = 1.0 / graph.size();
    node.value().c = 1.0 / graph.size();
  }
  
  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.001;
  double t_start = 0.0;
  double t_end   = 10.0;

  PointConstraint<GraphType> ptC = PointConstraint<GraphType>();
  PlaneConstraint<GraphType> pc = PlaneConstraint<GraphType>();
  // SphereConstraint<GraphType> sc = SphereConstraint<GraphType>();
  SphereBurnConstraint<GraphType> sbc = SphereBurnConstraint<GraphType>();
  auto allConstraint = make_combined_constraint(ptC, pc, sbc);

  for (double t = t_start; t < t_end; t += dt) {
    // symp_euler_step(graph, t, dt, Problem1Force());
    auto threeForce = make_combined_force(GravityForce(), MassSpringForce(), DampingForce());
    auto allForce = make_combined_force(threeForce, WindForce());
    symp_euler_step(graph, t, dt, allForce);
    // apply Point, Plane and Sphere constraints
    allConstraint(graph, t);

    viewer.clear();
    node_map.clear();
    // Update viewer with nodes' new positions
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.add_edges ( graph . edge_begin () , graph . edge_end () , node_map );

    viewer.set_label(t);

    // Make graphics more smoother :)
    //CS207::sleep(0.001);
  }

  return 0;
}
