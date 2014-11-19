#include "Point.hpp"
#include "Mesh.hpp"



typedef Mesh<int, bool, int> MeshType;
typedef MeshType::size_type size_type;


int main() {
    MeshType mesh;
    mesh.add_node(Point(-1,0,0));
    mesh.add_node(Point(0,0,0));
    mesh.add_node(Point(2,2,0));
    mesh.add_triangle(mesh.node(0), mesh.node(1), mesh.node(2));

    std::cout << "e_n1: " << mesh.triangle(0).edge(1).node1().position() << std::endl;
    std::cout << "e_n2: " << mesh.triangle(0).edge(1).node2().position() << std::endl;
    std::cout << "out norm of e: " << mesh.out_norm(mesh.triangle(0), mesh.triangle(0).edge(1)) << std::endl;
    std::cout << "area: " << mesh.triangle(0).area() << std::endl;
    return 0;
}