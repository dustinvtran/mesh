## Mesh

### Authors
* Ye Kuang \<yekuang@g.harvard.edu\>
* Dustin Tran \<dtran@g.harvard.edu\>

### Example
![](pond.gif)
```bash
make shallow_water
./shallow_water data/pond3.nodes data/pond3.tris 0
```
One can run the mesh on each of the three initial conditions by specifying `0`
(a pebble ripple), `1` (a sharp wave), or `2` (a dam break) for the third
argument of the `shallow_water` binary; by default it uses the `DamBreak`
condition.

### Description
We store a `Graph` object templated with node value type (`double` by default)
and edge value type `MeshEdgeValue`. `MeshEdgeValue` stores the left and right
triangles of an edge (chosen such that the left triangle corresponds to the one
left of the vector from node `i` to the node `j`, where `i<j`), and also the
template `E`, which is the edge value type one originally stores in a graph edge
(`double` by default).

We also store a vector of `internal_triangle` objects, whose index refers to a
triangle's index, and the internal triangle data refers to the associated
information for each triangle.

In order to iterate through all the triangles inside a `Mesh`, we define a class
called `TriangleIterator`. which stores a pointer `Mesh*` and a `index_`, which
points to the index of triangle being visited currently in that mesh.

As for the iteration of all the adjacent triangles of `Node` and `Triangle`, we
define a class called `NodeTriangleIterator` and `TriTriangleIterator`, which
stores a pointer `Mesh*` and one or two indices in order to keep track of the
iterative procedure during incrementing and dereferencing.

To avoid redundancy when setting up multiple initial conditions (and mostly for
fun), we implement a [Curiously Recursive Template
Pattern](http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern).
This constructs a base class from which each initial condition's class can
inherit from, allowing the initial conditions to merely state the changes done
onto node values, and then we run the same code to propagate this to the other
values.
