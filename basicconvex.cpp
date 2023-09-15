#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>  // for std::setw
#include <cmath>
#include <random>
#include <limits>
#include <unordered_set>
#include <stack>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstddef>
#include <new>
#include <memory>
#include <omp.h>
#include <deque>
#include <algorithm> // for std::all_of
#include <cstdlib> // for std::exit
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
// #include "/tinyobj_loader_opt.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Polyhedron_3<K> Polyhedron;
class Component;  // forward declaration

// Simplified Vector3D
struct Vector3D {
    double x, y, z;
    Vector3D() : x(0), y(0), z(0) {}
    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}
};

std::ostream& operator<<(std::ostream& os, const Vector3D& vec) {
    os << "(" << vec.x /* replace with your member variable */ << ", " << vec.y /* replace with your member variable */ << ", " << vec.z /* replace with your member variable */ << ")";
    return os;
}

class Triangle {
public:
    Vector3D vec1, vec2, vec3;
    Vector3D normal;
    double area;
    double e_min;
    double e_max;
    Vector3D color;

    Triangle() : vec1(), vec2(), vec3(), area(0.0), e_min(0.0), e_max(0.0) {}

    Triangle(const Vector3D& a, const Vector3D& b, const Vector3D& c);
};

Triangle::Triangle(const Vector3D& a, const Vector3D& b, const Vector3D& c)
    : vec1(a), vec2(b), vec3(c) {
    
}

// Overload the << operator for the Triangle class
std::ostream& operator<<(std::ostream& os, const Triangle& triangle) {
    os << "Triangle { ";
    os << "vec1: " << triangle.vec1 << ", ";
    os << "vec2: " << triangle.vec2 << ", ";
    os << "vec3: " << triangle.vec3 << ", ";
    os << "area: " << triangle.area << ", ";
    os << "e_min: " << triangle.e_min << ", ";
    os << "e_max: " << triangle.e_max;
    os << " }";
    return os;
}



class Hull {
public:
    std::vector<Point_3> points;
    Polyhedron tempP;

    Hull() {}
    Hull(const Component* component);


    static const long double EPSILON;  // Using long double for higher precision

    struct PointComparator {
        bool operator()(const Point_3& a, const Point_3& b) const {
            auto almostEqual = [](long double x, long double y, long double epsilon) {  // Changed to long double
                return std::abs(x - y) <= epsilon * std::max(1.0L, std::max(std::abs(x), std::abs(y)));  // Notice the "L" suffix for long double literals
            };
            if (!almostEqual((long double)a.x(), (long double)b.x(), EPSILON)) return a.x() < b.x();  // Cast to long double
            if (!almostEqual((long double)a.y(), (long double)b.y(), EPSILON)) return a.y() < b.y();
            if (!almostEqual((long double)a.z(), (long double)b.z(), EPSILON)) return a.z() < b.z();
            return false;  // Points are considered equal
        }
    };

    std::set<Point_3, PointComparator> pointSet;

    bool Hull::addUniquePoint(const Point_3& newPoint) {
        auto [iter, inserted] = pointSet.emplace(newPoint);
        if (inserted) {  
            points.push_back(newPoint);
        }
        return inserted;  
    }

    void Hull::add(const Triangle& triangle) {
        Point_3 p1(triangle.vec1.x, triangle.vec1.y, triangle.vec1.z);
        Point_3 p2(triangle.vec2.x, triangle.vec2.y, triangle.vec2.z);
        Point_3 p3(triangle.vec3.x, triangle.vec3.y, triangle.vec3.z);

        bool p1Added = addUniquePoint(p1);
        bool p2Added = addUniquePoint(p2);
        bool p3Added = addUniquePoint(p3);

        if (p1Added || p2Added || p3Added) { 
            computeHull();
        }
    }




    void computeHull() {
        if (points.size() < 4) {
            std::cout << "Not enough points for a valid 3D hull. Need at least 4." << std::endl;
            return;
        }
        
        std::cout << "Before computeHull, point count: " << points.size() << std::endl;
        
        // Clear the existing polyhedron
        CGAL::convex_hull_3(points.begin(), points.end(), tempP);

        std::cout << "Debug - Temp Hull Number of facets: " << tempP.size_of_facets() << std::endl;
        std::cout << "Debug - Temp Hull Number of vertices: " << tempP.size_of_vertices() << std::endl;

        // Debug: print the number of facets and vertices in the hull
        std::cout << "Number of facets in the hull: " << tempP.size_of_facets() << std::endl;
        std::cout << "Number of vertices in the hull: " << tempP.size_of_vertices() << std::endl;

        std::cout << "After computeHull, point count: " << points.size() << std::endl;
    }

    double tetrahedronVolume(const Point_3& A, const Point_3& B, const Point_3& C, const Point_3& D) {
        std::cout << "Calculating tetrahedron volume for points: " << A << ", " << B << ", " << C << ", " << D << std::endl;
        K::Vector_3 AB(B - A), AC(C - A), AD(D - A);

        // Calculate the cross product of AC and AD
        K::Vector_3 cross_product = CGAL::cross_product(AC, AD);

        // Compute the dot product of AB and the cross product of AC and AD
        double volume = std::abs(AB * cross_product) / 6.0;
        
        return volume;
    }



    
    double volume() const {
        if (tempP.is_empty()) {
            return 0.0;
        }

        double total_volume = 0.0;
        for (auto facet = tempP.facets_begin(); facet != tempP.facets_end(); ++facet) {
            auto h = facet->halfedge();
            Point_3 A = h->vertex()->point();
            Point_3 B = h->next()->vertex()->point();
            Point_3 C = h->opposite()->vertex()->point();
            Point_3 O(0, 0, 0); // Origin

            // Calculate tetrahedron volume
            double currentTetrahedronVolume = std::abs((A - O) * CGAL::cross_product(B - O, C - O)) / 6.0;

            std::cout << "Tetrahedron Volume: " << currentTetrahedronVolume << std::endl;
            total_volume += currentTetrahedronVolume;
        }

        return total_volume;
    }



    // Return the current set of points in the hull
    std::vector<Point_3> getPoints() const { 
        return points; 
    }

    // Return the current number of points in the hull
    int getPointCount() const { 
        return points.size(); 
    }

    void clear() {
        pointSet.clear();
        points.clear();
        tempP.clear();
    }
};

const long double Hull::EPSILON = 1e-9;  // Initialize outside the class definition

class Component {
public:
    std::vector<Triangle> triangles;
    Hull hull;

    Component() : hull() {}

    void addTriangleAndUpdateHull(const Triangle& triangle) {
        triangles.push_back(triangle);
        hull.add(triangle);  // No need to call computeHull here, as it's already done in add.
    }

    // Function to choose the next triangle by evaluating the hull
    std::pair<Triangle, double> chooseNextTriangle(const std::vector<Triangle>& remainingTriangles) {
        Triangle bestTriangle;
        double bestScore = std::numeric_limits<double>::max();

        for (const Triangle& triangle : remainingTriangles) {
            Hull tempHull = hull; // Making a copy of the existing hull

            tempHull.add(triangle);  // add will internally manage whether to recompute the hull or not

            if (tempHull.getPointCount() >= 4) {
                double newVolume = tempHull.volume();
                double oldVolume = hull.volume();
                double volumeChange = std::abs(newVolume - oldVolume);

                if (volumeChange < bestScore) {
                    bestScore = volumeChange;
                    bestTriangle = triangle;
                }
            }
        }

        return {bestTriangle, bestScore};
    }
};


Hull::Hull(const Component* component) {
    if (component) {
        for (const auto& triangle : component->triangles) {
            add(triangle);
        }
        computeHull();
    }
}


struct TriangleIndices {
    int v1, v2, v3;

    // Default constructor
    TriangleIndices() : v1(0), v2(0), v3(0) {}

    // Parameterized constructor
    TriangleIndices(int v1, int v2, int v3) : v1(v1), v2(v2), v3(v3) {}
};

std::vector<Triangle> indicesToTriangles(const std::vector<TriangleIndices>& indicesList, const std::vector<Vector3D>& vertices) {
    std::vector<Triangle> triangles;
    for (const auto& indices : indicesList) {
        Vector3D v1 = vertices[indices.v1];
        Vector3D v2 = vertices[indices.v2];
        Vector3D v3 = vertices[indices.v3];

        Triangle tri(v1, v2, v3);
        triangles.push_back(tri);
    }
    return triangles;
}

std::pair<Triangle, std::vector<double>> chooseNextTriangleByHull(
    const Component& tempComponent, 
    const std::vector<Triangle>& remainingTriangles
) {
    Triangle bestTriangle;
    double bestScore = std::numeric_limits<double>::max();
    std::vector<double> bestVolumeChanges;

    std::cout << "Initial hull volume: " << tempComponent.hull.volume() << std::endl;

    for (const Triangle& triangle : remainingTriangles) {
            Hull tempHull = tempComponent.hull; // Making a copy of the hull

            tempHull.add(triangle);  // add will internally manage whether to recompute the hull or not

            if (tempHull.getPointCount() >= 4) {
            double newVolume = tempHull.volume();
            std::cout << "Volume immediately after adding triangle: " << newVolume << std::endl;
            double oldVolume = tempComponent.hull.volume();
            double volumeChange = std::abs(newVolume - oldVolume);

            std::cout << "Old Volume: " << oldVolume << ", New Volume: " << newVolume 
                      << ", Volume Change: " << volumeChange << std::endl;

            if (volumeChange < bestScore) {
                bestScore = volumeChange;
                bestTriangle = triangle;
                bestVolumeChanges.clear();
                bestVolumeChanges.push_back(volumeChange);
            }
        }
    }

    std::cout << "Best Score: " << bestScore << ", Best Triangle: " << bestTriangle << std::endl;
    return {bestTriangle, bestVolumeChanges};
}


int main() {
    std::cout << "Starting program..." << std::endl;
    // std::vector<Triangle> allTriangles = {};  // Fill with triangles from your data source
    std::vector<double> weights = {.3, .3, .3, 0};
    std::vector<Vector3D> allVertices;    
    double tau_S =.01;  // Threshold for shape similarity
    std::string inputfile = "/home/kingpin/Documents/blender/cottage.obj";
    std::string err;
    // Read the OBJ file into a buffer
    std::ifstream objFile(inputfile, std::ios::in | std::ios::binary | std::ios::ate);
    std::streamsize size = objFile.tellg();
    objFile.seekg(0, std::ios::beg);

    std::vector<char> objBuffer(size);
    if (!objFile.is_open()) {
        std::cerr << "Failed to open the OBJ file." << std::endl;
        return -1;
    }
    if (!objFile.read(objBuffer.data(), size)) {
        std::cerr << "Failed to read the OBJ file." << std::endl;
        return -1;
    }

    std::cout << "Finished reading OBJ file." << std::endl;


    // Parse the OBJ
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    std::string warn;

    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str());

    if (!warn.empty()) {
        std::cout << "WARN: " << warn << std::endl;
    }

    if (!err.empty()) {
        std::cerr << err << std::endl;
    }

    if (!ret) {
        std::cerr << "Failed to load/parse .obj file!" << std::endl;
        return -1;
    }

    std::cout << "Total vertices in OBJ: " << attrib.vertices.size()/3 << std::endl;

    for (size_t i = 0; i < attrib.vertices.size(); i += 3) {
        Vector3D v(attrib.vertices[i], attrib.vertices[i + 1], attrib.vertices[i + 2]);
        // std::cout << "Reading Vertex: (" << attrib.vertices[i] << "," << attrib.vertices[i + 1] << "," << attrib.vertices[i + 2] << ")" << std::endl;
        allVertices.push_back(v);
    }

    // Instead of converting vertices to Triangle structures, we will simply store their indices for writing back.
    std::vector<TriangleIndices> allTriangleIndices;

    std::cout << "Number of shapes: " << shapes.size() << std::endl;

    for (const auto& shape : shapes) {
        const auto& mesh = shape.mesh;

        size_t index_offset = 0;
            for (size_t f = 0; f < mesh.num_face_vertices.size(); f++) {
                int fv = mesh.num_face_vertices[f];

                if (fv == 3) {
                    TriangleIndices tri;
                    tri.v1 = mesh.indices[index_offset + 0].vertex_index;
                    tri.v2 = mesh.indices[index_offset + 1].vertex_index;
                    tri.v3 = mesh.indices[index_offset + 2].vertex_index;
                    allTriangleIndices.push_back(tri);
                }
                else if (fv == 4) {
                    // Handle quads by splitting them into two triangles
                    TriangleIndices tri1;
                    tri1.v1 = mesh.indices[index_offset + 0].vertex_index;
                    tri1.v2 = mesh.indices[index_offset + 1].vertex_index;
                    tri1.v3 = mesh.indices[index_offset + 2].vertex_index;
                    allTriangleIndices.push_back(tri1);

                    TriangleIndices tri2;
                    tri2.v1 = mesh.indices[index_offset + 0].vertex_index;
                    tri2.v2 = mesh.indices[index_offset + 2].vertex_index;
                    tri2.v3 = mesh.indices[index_offset + 3].vertex_index;
                    allTriangleIndices.push_back(tri2);
                }
                else {
                    std::cerr << "Warning: face with " << fv << " vertices. Not handled in this code." << std::endl;
                }
                
                index_offset += fv;
            }
    }

    std::cout << "First 10 triangle indices:" << std::endl;
    for (int i = 0; i < 10 && i < allTriangleIndices.size(); i++) {
        std::cout << "Triangle Indices " << i << ": " << allTriangleIndices[i].v1 << ", " << allTriangleIndices[i].v2 << ", " << allTriangleIndices[i].v3 << std::endl;
    }
    
    std::vector<Triangle> allTriangles = indicesToTriangles(allTriangleIndices, allVertices);

    std::map<std::string, Component> component_hulls;

    // Initialize a single component as an example
    Component initial_component;
    initial_component.addTriangleAndUpdateHull(allTriangles[0]);
    component_hulls["initial_component"] = initial_component;

    std::vector<std::vector<Triangle>> candidate_triangle_sets;  
    // Assuming you populate this vector with your sets of candidate triangles.

    // Add your new logic here
    for (auto& triangle_set : candidate_triangle_sets) {
        std::map<std::string, Triangle> best_triangle_per_component;

        // Loop through each component
        for (auto& [component_name, component] : component_hulls) {
            auto [best_triangle, best_score] = component.chooseNextTriangle(triangle_set);
            best_triangle_per_component[component_name] = best_triangle;
        }

        // Add best triangles to each component's hull
        for (auto& [component_name, best_triangle] : best_triangle_per_component) {
            component_hulls[component_name].addTriangleAndUpdateHull(best_triangle);
        }
    }

    // Test 3: Check points and volume after adding a second triangle
    // std::cout << "Hull Volume after adding 2 triangles: " << component.hull.volume() << std::endl;

    return 0;
}
