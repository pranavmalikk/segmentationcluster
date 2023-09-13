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

// Simplified Vector3D
struct Vector3D {
    double x, y, z;
    Vector3D() : x(0), y(0), z(0) {}
    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}
};

// Simplified Triangle
struct Triangle {
    Vector3D v1, v2, v3;
    Triangle(Vector3D v1, Vector3D v2, Vector3D v3) : v1(v1), v2(v2), v3(v3) {}
};

// Simplified Component
struct Component {
    std::vector<Triangle> triangles;
    // Hull is simplified to a single point for demonstration
    Vector3D hullPoint;

    // Updates the "hull" when a new triangle is added
    void updateHull() {
        // For demonstration, let's say the hull point is just the average of all triangle vertices
        Vector3D sum;
        for(const auto& t : triangles) {
            sum.x += (t.v1.x + t.v2.x + t.v3.x) / 3.0;
            sum.y += (t.v1.y + t.v2.y + t.v3.y) / 3.0;
            sum.z += (t.v1.z + t.v2.z + t.v3.z) / 3.0;
        }
        hullPoint = {sum.x / triangles.size(), sum.y / triangles.size(), sum.z / triangles.size()};
    }
};

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

// Function to choose the next triangle based on the hull (simplified)
Triangle chooseNextTriangleByHull(const Component& component, const std::vector<Triangle>& candidates) {
    // Just select the triangle whose first vertex is closest to the hull point
    double minDist = std::numeric_limits<double>::max();
    Triangle bestTriangle = candidates[0];
    for(const auto& t : candidates) {
        double dist = std::sqrt(std::pow(t.v1.x - component.hullPoint.x, 2) +
                                std::pow(t.v1.y - component.hullPoint.y, 2) +
                                std::pow(t.v1.z - component.hullPoint.z, 2));
        if(dist < minDist) {
            minDist = dist;
            bestTriangle = t;
        }
    }
    return bestTriangle;
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

    // Initialize a component with one triangle
    Component component;
    component.triangles.push_back(allTriangles[0]);
    component.updateHull();

    // Choose the next triangle to add to the component
    Triangle nextTriangle = chooseNextTriangleByHull(component, allTriangles);

    // Add the chosen triangle to the component and update the hull
    component.triangles.push_back(nextTriangle);
    component.updateHull();

    // Output for verification
    std::cout << "New Hull Point: (" << component.hullPoint.x << ", " << component.hullPoint.y << ", " << component.hullPoint.z << ")\n";

    return 0;
}
