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
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/convex_hull_3_to_face_graph.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Surface_mesh.h>
#include <list>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstddef>
#include <new>
#include <memory>
#include <omp.h>
#include <deque>
#include <queue> // for std::priority_queue
#include <utility> // for std::pair
#include <algorithm> // for std::all_of
#include <cstdlib> // for std::exit
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
// #include "/tinyobj_loader_opt.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
class Triangle;
class TriangleAdjacency;
class Hull; 
// Assuming a Vector3D structure is present somewhere, used for normals.
class Vector3D {
public:
    double x, y, z;

    // Constructors
    Vector3D() : x(0), y(0), z(0) {}
    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}

    // Subtraction
    Vector3D operator-(const Vector3D& other) const {
        return Vector3D(x - other.x, y - other.y, z - other.z);
    }

    // Cross product
    Vector3D cross(const Vector3D& other) const {
        double newX = y * other.z - z * other.y;
        double newY = z * other.x - x * other.z;
        double newZ = x * other.y - y * other.x;
        return Vector3D(newX, newY, newZ);
    }

    Vector3D operator/(double scalar) const {
        if (std::abs(scalar) < 1e-9) { // Avoid division by zero
            return *this; // Return the original vector or handle this case differently
        }
        return Vector3D(x / scalar, y / scalar, z / scalar);
    }

    // Overload the += operator
    Vector3D& operator+=(const Vector3D& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this; // Return this instance
    }

    // Overload the /= operator for scalar division
    Vector3D& operator/=(double scalar) {
        if (scalar != 0.0) {  // Avoid division by zero
            x /= scalar;
            y /= scalar;
            z /= scalar;
        }
        return *this; // Return this instance
    }

    bool operator<(const Vector3D& rhs) const {
        if (x != rhs.x) return x < rhs.x;
        if (y != rhs.y) return y < rhs.y;
        return z < rhs.z;
    }

    Vector3D operator+(const Vector3D& other) const {
        Vector3D result;
        result.x = this->x + other.x;
        result.y = this->y + other.y;
        result.z = this->z + other.z;
        return result;
    }

    double dot(const Vector3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    // Length (magnitude) of the vector
    double length() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    // Check if vector is a zero vector
    bool isZero() const {
        const double epsilon = 1e-9; // Adjust this to your needs
        return std::abs(x) < epsilon && std::abs(y) < epsilon && std::abs(z) < epsilon;
    }

    bool operator==(const Vector3D& other) const {
        return x == other.x && y == other.y && z == other.z;
    }

    std::string toString() const {
        return "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")";
    }

    void normalize() {
        double len = length();
        if (len > 0) {
            x /= len;
            y /= len;
            z /= len;
        }
    }
};

struct Vector3DHash {
    std::size_t operator()(const Vector3D& vec) const {
        std::size_t hx = std::hash<double>()(vec.x);
        std::size_t hy = std::hash<double>()(vec.y);
        std::size_t hz = std::hash<double>()(vec.z);
        return hx ^ (hy << 1) ^ hz;
    }
};

namespace std {
    template <>
    struct hash<Vector3D> {
        std::size_t operator()(const Vector3D& vec) const {
            Vector3DHash hasher;
            return hasher(vec);
        }
    };
}

class Vertex {
public:
    int x, y, z;  // 3D coordinates

    Vertex() : x(0), y(0), z(0) {}
    Vertex(int index) : x(index), y(index), z(index) {}  // New constructor
    Vertex(int x, int y, int z) : x(x), y(y), z(z) {}
    Vertex(const Vector3D& vec) : x(vec.x), y(vec.y), z(vec.z) {}

    bool operator==(const Vertex& other) const {
        return x == other.x && y == other.y && z == other.z;
    }

    bool operator!=(const Vertex& other) const {
        return !(*this == other);
    }

    bool isValid() const {
        return !std::isnan(x) && !std::isnan(y) && !std::isnan(z);
    }

    Vector3D toVector3D() const {
        return Vector3D(x, y, z);  // assuming Vertex has x, y, and z as public members.
    }

    std::string toString() const {
        return "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")";
    }
};

namespace std {
    template <> struct hash<Vertex> {
        size_t operator()(const Vertex& v) const {
            return hash<int>()(v.x) ^ hash<int>()(v.y) ^ hash<int>()(v.z);
        }
    };
}

struct TriangleEdge {
    Vector3D v1, v2;

    TriangleEdge(const Vector3D& a, const Vector3D& b) : v1(a), v2(b) {
        if (v2 < v1) {
            std::swap(v1, v2);
        }
    }

    bool operator==(const TriangleEdge& other) const {
        return v1 == other.v1 && v2 == other.v2;
    }

    bool operator!=(const TriangleEdge& other) const {
        return !(*this == other);
    }

    bool operator<(const TriangleEdge& other) const {
        if (v1 == other.v1) return v2 < other.v2;
        return v1 < other.v1;
    }

    std::string toString() const {
        return "[" + v1.toString() + ", " + v2.toString() + "]";
    }
};

struct EdgeHash {
    Vector3DHash vecHash;  // Using your custom hash function for Vector3D

    size_t operator()(const TriangleEdge& edge) const {
        size_t h1 = vecHash(edge.v1);
        size_t h2 = vecHash(edge.v2);
        return h1 ^ (h2 << 1); // XOR and shift h2 before combining with h1 for better distribution
    }
};

class Triangle {
public:
    Vector3D vec1, vec2, vec3;
    Vector3D normal;
    double area;
    double e_min;
    double e_max;
    Vector3D color;
    std::string type; // The type of this triangle

    Triangle() : area(0), e_min(0), e_max(0) {}

    Triangle(const Vector3D& a, const Vector3D& b, const Vector3D& c);

    bool isValid() const;
    bool operator==(const Triangle& other) const;
    bool operator!=(const Triangle& other) const;
    bool operator<(const Triangle& other) const;
    void setColor(const Vector3D& assignedColor);
    std::vector<TriangleEdge> getEdges() const;
    bool isAdjacent(const Triangle& other, const class TriangleAdjacency& adjacency) const;
    bool isNewAdjacent(const Triangle& other, const class TriangleAdjacency& adjacency) const;  // Added this
    std::string toString() const;  // Added toString() method
    bool isDegenerate() const;

};

namespace std {
    template <>
    struct hash<Triangle> {
        size_t operator()(const Triangle& t) const {
            std::hash<Vector3D> vertexHasher;
            size_t h1 = vertexHasher(t.vec1);
            size_t h2 = vertexHasher(t.vec2);
            size_t h3 = vertexHasher(t.vec3);
            return h1 ^ (h2 << 1) ^ (h3 << 2);  // or use boost::hash_combine (better)
        }
    };
}

class TriangleAdjacency {
private:
    typedef TriangleEdge Edge;
    std::unordered_map<Edge, std::unordered_set<Triangle>, EdgeHash> adjacencyMap;


public:
    void addTriangle(const Triangle& triangle);
    bool isAdjacent(const Triangle& t1, const Triangle& t2) const;
    bool isNewAdjacent(const Triangle& t1, const Triangle& t2) const;
    std::string toString() const;  // Added toString() method
    std::unordered_set<Triangle> getAdjacentTriangles(const Triangle& t) const; 
};

// Triangle methods
Triangle::Triangle(const Vector3D& a, const Vector3D& b, const Vector3D& c)
    : vec1(a), vec2(b), vec3(c) {
    Vector3D edge1 = vec2 - vec1;
    Vector3D edge2 = vec3 - vec1;
    normal = edge1.cross(edge2);
    double normalLength = normal.length();
    normal = normal / normalLength;
    if (normalLength > 1e-9) {
        normal = normal / normalLength;  // This line normalizes the normal vector
    }
    area = 0.5 * normalLength;
    double len1 = edge1.length();
    double len2 = edge2.length();
    double len3 = (vec3 - vec2).length();
    e_min = std::min({len1, len2, len3});
    e_max = std::max({len1, len2, len3});
    
}

bool Triangle::isValid() const {
    return !vec1.isZero() || !vec2.isZero() || !vec3.isZero();
}

bool Triangle::operator==(const Triangle& other) const {
    return (vec1 == other.vec1 && vec2 == other.vec2 && vec3 == other.vec3) ||
           (vec1 == other.vec1 && vec2 == other.vec3 && vec3 == other.vec2) ||
           (vec1 == other.vec2 && vec2 == other.vec1 && vec3 == other.vec3) ||
           (vec1 == other.vec2 && vec2 == other.vec3 && vec3 == other.vec1) ||
           (vec1 == other.vec3 && vec2 == other.vec1 && vec3 == other.vec2) ||
           (vec1 == other.vec3 && vec2 == other.vec2 && vec3 == other.vec1);
}


bool Triangle::operator!=(const Triangle& other) const {
    return !(*this == other);
}

bool Triangle::operator<(const Triangle& other) const {
    return this->area < other.area;
}

void Triangle::setColor(const Vector3D& assignedColor) {
    color = assignedColor;
}

std::vector<TriangleEdge> Triangle::getEdges() const {
    return {TriangleEdge(vec1, vec2), TriangleEdge(vec2, vec3), TriangleEdge(vec3, vec1)};
}

bool Triangle::isAdjacent(const Triangle& other, const TriangleAdjacency& adjacency) const {
    return adjacency.isAdjacent(*this, other);
}

bool Triangle::isNewAdjacent(const Triangle& other, const TriangleAdjacency& adjacency) const {
    return adjacency.isNewAdjacent(*this, other);
}

std::string Triangle::toString() const {
    return "Triangle: [" + vec1.toString() + ", " + vec2.toString() + ", " + vec3.toString() + "]";
}

bool Triangle::isDegenerate() const {
    double area = 0.5 * std::abs((vec1.x - vec2.x)*(vec1.y - vec3.y) - (vec1.x - vec3.x)*(vec1.y - vec2.y));
    return area == 0;
}


// TriangleAdjacency methods
void TriangleAdjacency::addTriangle(const Triangle& triangle) {
    std::vector<TriangleEdge> edges = {
        TriangleEdge(triangle.vec1, triangle.vec2),
        TriangleEdge(triangle.vec2, triangle.vec3),
        TriangleEdge(triangle.vec3, triangle.vec1)
    };

    for (const TriangleEdge& edge : edges) {
        adjacencyMap[edge].insert(triangle);
    }
}

bool TriangleAdjacency::isAdjacent(const Triangle& t1, const Triangle& t2) const {
    std::vector<TriangleEdge> edges = {
        TriangleEdge(t1.vec1, t1.vec2),
        TriangleEdge(t1.vec2, t1.vec3),
        TriangleEdge(t1.vec3, t1.vec1)
    };

    for (const TriangleEdge& edge : edges) {
        if (adjacencyMap.count(edge) && adjacencyMap.at(edge).count(t2)) {
            return true;
        }
    }
    return false;
}

bool TriangleAdjacency::isNewAdjacent(const Triangle& t1, const Triangle& t2) const {
    // Implement your new adjacency logic here
    
    // Example: Share at least two vertices
    std::unordered_set<Vector3D> vertices1 = {t1.vec1, t1.vec2, t1.vec3};
    std::unordered_set<Vector3D> vertices2 = {t2.vec1, t2.vec2, t2.vec3};
    std::vector<Vector3D> commonVertices;
    
    for (const auto& vertex : vertices1) {
        if (vertices2.count(vertex)) {
            commonVertices.push_back(vertex);
        }
    }
    
    if (commonVertices.size() >= 2) {
        return true;
    }

    return false; // Replace this with your logic
}

std::string TriangleAdjacency::toString() const {
    std::ostringstream oss;
    oss << "TriangleAdjacency: {" << std::endl;
    
    for (const auto& entry : adjacencyMap) {
        const Edge& edge = entry.first;
        const std::unordered_set<Triangle>& adjacentTriangles = entry.second;
        
        oss << "  Edge: " << edge.toString() << " -> Adjacent Triangles: {";
        
        for (const Triangle& triangle : adjacentTriangles) {
            oss << triangle.toString() << ", ";
        }
        
        oss << "}," << std::endl;
    }
    
    oss << "}" << std::endl;
    return oss.str();
}

std::unordered_set<Triangle> TriangleAdjacency::getAdjacentTriangles(const Triangle& t) const {
    std::unordered_set<Triangle> adjacentTriangles;
    for (const TriangleEdge& edge : t.getEdges()) {
        auto it = adjacencyMap.find(edge);
        if (it != adjacencyMap.end()) {
            adjacentTriangles.insert(it->second.begin(), it->second.end());
        }
    }
    adjacentTriangles.erase(t);  // Remove the triangle itself from its adjacent list
    return adjacentTriangles;
}


struct TriangleIndices {
    int v1, v2, v3;

    // Default constructor
    TriangleIndices() : v1(0), v2(0), v3(0) {}

    // Parameterized constructor
    TriangleIndices(int v1, int v2, int v3) : v1(v1), v2(v2), v3(v3) {}
};

Vector3D computeNormal(const Vector3D& v1, const Vector3D& v2, const Vector3D& v3) {
    Vector3D edge1 = v2 - v1;
    Vector3D edge2 = v3 - v1;
    Vector3D normal = edge1.cross(edge2);
    normal.normalize(); // Assuming you have a normalize method
    return normal;
}

std::vector<Triangle> indicesToTriangles(const std::vector<TriangleIndices>& indices, const std::vector<Vector3D>& vertices) {
    std::vector<Triangle> triangles;
    for (const auto& index : indices) {
        Vector3D v1 = vertices[index.v1];
        Vector3D v2 = vertices[index.v2];
        Vector3D v3 = vertices[index.v3];
        Triangle t(v1, v2, v3);
        
        // Compute the normal and normalize it
        t.normal = computeNormal(v1, v2, v3);

        // Compute the magnitude of the normal vector manually
        double magnitude = std::sqrt(t.normal.x * t.normal.x + t.normal.y * t.normal.y + t.normal.z * t.normal.z);

        // Check for a degenerate triangle based on the magnitude of the normal vector
        if (magnitude < 1e-6) {  // Adjust threshold as needed
            std::cerr << "Degenerate triangle detected, skipping." << std::endl;
            continue;  // Skip this triangle and continue with the next iteration
        }

        triangles.push_back(t);
    }
    return triangles;
}


struct Cluster {
    std::vector<Triangle> triangles;

    bool operator==(const Cluster& other) const {
        return triangles == other.triangles;
    }
};

struct TupleHash {
    std::size_t operator()(const std::tuple<int, int, int>& k) const {
    size_t h1 = std::hash<int>{}(std::get<0>(k));
    size_t h2 = std::hash<int>{}(std::get<1>(k));
    size_t h3 = std::hash<int>{}(std::get<2>(k));
    return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

struct IndexedTriangle {
    int clusterIndex;
    Triangle triangle;
};

class SpatialHash {
private:
    double cellSize;
    double invCellSize;  // Reciprocal of cellSize for optimization
    std::unordered_map<std::tuple<int, int, int>, std::vector<IndexedTriangle>, TupleHash> hashTable;
    std::unordered_map<Triangle, std::vector<Triangle>> adjacencyMap; 
    std::unordered_map<Triangle, std::array<std::tuple<int, int, int>, 3>> precomputedHashes;
    // std::unordered_map<Triangle, std::vector<Triangle>> precomputedNeighbors;

    std::tuple<int, int, int> hash(const Vector3D& vec) {
        int x = static_cast<int>(std::floor(vec.x * invCellSize));
        int y = static_cast<int>(std::floor(vec.y * invCellSize));
        int z = static_cast<int>(std::floor(vec.z * invCellSize));

        return {x, y, z};
    }
    

    std::array<std::tuple<int, int, int>, 3> getTriangleHashes(const Triangle& t) {
        return {hash(t.vec1), hash(t.vec2), hash(t.vec3)};
    }

public:
    SpatialHash(double size) : cellSize(size), invCellSize(1.0 / size) {}


    std::vector<std::string> keys() const {
        std::vector<std::string> allKeys;
        for (const auto& [key, _] : hashTable) {
            auto [x, y, z] = key;
            allKeys.push_back(std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z));
        }
        return allKeys;
    }

    size_t size() const {
        return hashTable.size();
    }

    std::unordered_set<int> getNeighboringClusters(const Triangle& t, const TriangleAdjacency& adjacency) {
        std::unordered_set<int> clusterIndices;

        // Use the iterator from find to avoid a second lookup
        auto it = precomputedHashes.find(t);
        if (it == precomputedHashes.end()) {
            std::cerr << "Triangle not found in precomputed hashes: " << t.toString() << std::endl;
            return clusterIndices;
        }

        const auto& hashes = it->second;
        for (const auto& h : hashes) {
            // Use iterator for the hashTable lookup
            auto hashIt = hashTable.find(h);
            if (hashIt == hashTable.end()) continue;

            for (const auto& indexedTriangle : hashIt->second) {
                const Triangle& potentialNeighbor = indexedTriangle.triangle;
                if (potentialNeighbor != t) {
                    // Add an adjacency check here
                    // if (t.isAdjacent(potentialNeighbor, adjacency)) {
                        clusterIndices.insert(indexedTriangle.clusterIndex);
                    // }
                }
            }
        }
        return clusterIndices;
    }


    std::vector<Triangle> getPotentialNeighbors(const Triangle& t, const TriangleAdjacency& adjacency) {
        std::unordered_set<Triangle> neighbors;
        auto it = precomputedHashes.find(t);
        
        if (it == precomputedHashes.end()) {
            std::cout << "Error: Triangle hashes not precomputed for: " << t.toString() << std::endl;
            return {};
        }

        const auto& hashes = it->second;
        for (const auto& h : hashes) {
            for (const auto& indexedTriangle : hashTable[h]) {
                const Triangle& potentialNeighbor = indexedTriangle.triangle;
                if (potentialNeighbor != t && t.isAdjacent(potentialNeighbor, adjacency)) {
                    neighbors.insert(potentialNeighbor);
                }
            }
        }
        return std::vector<Triangle>(neighbors.begin(), neighbors.end());
    }

    std::unordered_set<int> getNeighboringClustersForCluster(const Cluster& cluster, const TriangleAdjacency& adjacency) {
        std::unordered_set<int> clusterIndices;
        std::unordered_set<Triangle> processedTriangles;
            
        for (const auto& triangle : cluster.triangles) {
            auto neighboringClustersForTriangle = getNeighboringClusters(triangle, adjacency);
            clusterIndices.insert(neighboringClustersForTriangle.begin(), neighboringClustersForTriangle.end());
        }
            
        return clusterIndices;
    }

    void precomputeTriangleHashes(const std::vector<Cluster>& clusters) {
        size_t totalTriangles = 0;
        for (size_t clusterIndex = 0; clusterIndex < clusters.size(); ++clusterIndex) {
            const Cluster& cluster = clusters[clusterIndex];
            for (const Triangle& triangle : cluster.triangles) {
                auto hashes = getTriangleHashes(triangle);
                precomputedHashes[triangle] = hashes;
                totalTriangles++;

                for (const auto& h : hashes) {
                    hashTable[h].emplace_back(IndexedTriangle{static_cast<int>(clusterIndex), triangle});
                }
            }
        }

    }


    void clear() {
        hashTable.clear();
    }
};

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

// Checks if two triangles share a vertex
bool shareVertex(const Triangle& a, const Triangle& b) {
    return a.vec1 == b.vec1 || a.vec1 == b.vec2 || a.vec1 == b.vec3 ||
           a.vec2 == b.vec1 || a.vec2 == b.vec2 || a.vec2 == b.vec3 ||
           a.vec3 == b.vec1 || a.vec3 == b.vec2 || a.vec3 == b.vec3;
}

// Depth First Search function
void DFS(int u, std::vector<std::vector<int>>& adjList, std::vector<bool>& visited) {
    // std::cout << "[Debug] Visiting vertex: " << u << std::endl;
    visited[u] = true;
    for (int v : adjList[u]) {
        if (!visited[v]) {
            DFS(v, adjList, visited);
        }
    }
}

bool areTrianglesContiguous(const std::vector<Triangle>& triangles) {
    int n = triangles.size();
    if (n == 0) return true;

    std::vector<std::vector<int>> adjList(n);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (shareVertex(triangles[i], triangles[j])) {
                adjList[i].push_back(j);
                adjList[j].push_back(i);
            }
        }
    }

    std::vector<bool> visited(n, false);
    int connectedComponents = 0;

    for (int i = 0; i < n; ++i) {
        if (!visited[i]) {
            // std::cout << "[Debug] Running DFS starting from: " << i << std::endl;
            DFS(i, adjList, visited);
            connectedComponents++;
        }
    }

    // std::cout << "[Debug] Number of connected components: " << connectedComponents << std::endl;

    return connectedComponents == 1;
}

void print_convex_hull(const CGAL::Surface_mesh<Point_3>& P) {
    std::cout << "Vertices of the convex hull:" << std::endl;
    for (auto v : P.vertices()) {
        Point_3 p = P.point(v);
        std::cout << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")" << std::endl;
    }

    std::cout << "Edges of the convex hull:" << std::endl;
    for (auto e : P.edges()) {
        auto v1 = P.vertex(e, 0);
        auto v2 = P.vertex(e, 1);
        Point_3 p1 = P.point(v1);
        Point_3 p2 = P.point(v2);
        std::cout << "[(" << p1.x() << ", " << p1.y() << ", " << p1.z() << ") - ("
                  << p2.x() << ", " << p2.y() << ", " << p2.z() << ")]" << std::endl;
    }

    // std::cout << "Faces of the convex hull:" << std::endl;
    // for (auto f : P.faces()) {
    //     std::cout << "Face: ";
    //     auto h = P.halfedge(f);  // Obtain a halfedge handle from face f
    //     if(h == CGAL::Surface_mesh<Point_3>::null_halfedge()) {
    //         std::cerr << "Error: null halfedge encountered." << std::endl;
    //         continue;  // Skip this face if we got a null halfedge
    //     }
    //     CGAL::Halfedge_around_face_circulator<CGAL::Surface_mesh<Point_3>> hf_circ(h, P), hf_end;
    //     do {
    //         auto v = P.target(*hf_circ);
    //         Point_3 p = P.point(v);
    //         std::cout << "(" << p.x() << ", " << p.y() << ", " << p.z() << ") ";
    //     } while (++hf_circ != hf_end);
    //     std::cout << std::endl;
    // }
}

void print_2d_convex_hull(const Polygon_2& P_2D) {
    std::cout << "Vertices of the 2D convex hull:" << std::endl;
    for (auto vertex_it = P_2D.vertices_begin(); vertex_it != P_2D.vertices_end(); ++vertex_it) {
        CGAL::Point_2<K> p = *vertex_it;
        std::cout << "(" << p.x() << ", " << p.y() << ")" << std::endl;
    }

    std::cout << "Edges of the 2D convex hull:" << std::endl;
    for (auto edge_it = P_2D.edges_begin(); edge_it != P_2D.edges_end(); ++edge_it) {
        CGAL::Segment_2<K> segment = *edge_it;
        CGAL::Point_2<K> source = segment.source();
        CGAL::Point_2<K> target = segment.target();
        std::cout << "[(" << source.x() << ", " << source.y() << ") - ("
                  << target.x() << ", " << target.y() << ")]" << std::endl;
    }
}


class Hull {
private:
    Delaunay T;  // Delaunay triangulation
    std::vector<Point_3> points;
    Surface_mesh P;  // Adjusted to Surface_mesh
    mutable double cachedVolume = -1;
    mutable bool isVolumeCacheValid = false;

public:

    Hull() = default;  // Add a default constructor if it's missing

    Hull(const std::vector<Triangle>& triangles) {
        for (const Triangle& triangle : triangles) {
            add(triangle);
        }
        computeHull();
    }

    // Copy constructor
    Hull(const Hull& other) {
        points = other.points;
        P = other.P;
    }

    Hull(const Hull& h1, const Hull& h2) {
        points = h1.points;
        points.insert(points.end(), h2.points.begin(), h2.points.end());
        computeHull();
    }

    // Add this method to access the private member P
    const Surface_mesh& getSurfaceMesh() const {
        return P;
    }

    std::string toString() const {
        std::ostringstream oss;
        // oss << "Hull with " << points.size() << " points and ";
        if (P.is_empty()) {
            oss << "no surface mesh constructed.";
        } else {
            oss << "surface mesh with " << P.number_of_vertices() << " vertices, "
                << P.number_of_halfedges() << " halfedges, and "
                << P.number_of_faces() << " facets.";
        }
        return oss.str();
    }

    bool arePointsCoplanar(const std::vector<Point_3>& points) {
        if (points.size() < 4) return true;  // Less than 4 points are always coplanar
        const Point_3& p1 = points[0];
        const Point_3& p2 = points[1];
        const Point_3& p3 = points[2];
        CGAL::Vector_3 normal = CGAL::cross_product(p2 - p1, p3 - p1);
        for (size_t i = 3; i < points.size(); ++i) {
            CGAL::Vector_3 vec = points[i] - p1;
            if (CGAL::scalar_product(normal, vec) != 0) {
                return false;
            }
        }
        return true;
    }

    void add(const Triangle& triangle) {
        addPoint(triangle.vec1);
        addPoint(triangle.vec2);
        addPoint(triangle.vec3);
    }

    void addPoint(const Vector3D& vec) {
        Point_3 point(vec.x, vec.y, vec.z);
        // T.insert(point);  // Insert point into Delaunay triangulation
        points.push_back(point);  // Store point in vector
    }

    void removePoint(const Point_3& point) {
        auto vertex_handle = T.nearest_vertex(point);
        if (vertex_handle != nullptr && vertex_handle->point() == point) {
            T.remove(vertex_handle);
        }
    }

    void computeHull() {
        try {
            CGAL::convex_hull_3(points.begin(), points.end(), P);
        } catch (const std::exception& e) {
            std::cerr << "Exception: " << e.what() << std::endl;
            exit(1);
        }
    }

    bool arePointsUnique() const {
        std::set<Point_3> uniquePoints(points.begin(), points.end());
        return uniquePoints.size() == points.size();
    }

    bool isEmpty() const {
        return P.is_empty();
    }

    bool isClosed() const {
        return CGAL::is_closed(P);
    }

    void updateVolume() const {  // add const qualifier here
        if (P.is_empty() || !CGAL::is_closed(P)) {
            cachedVolume = 0.0;
        } else {
            cachedVolume = CGAL::Polygon_mesh_processing::volume(P);
        }
        isVolumeCacheValid = true;  // set the cache flag to valid
    }

    double volume() const {
        if (!isVolumeCacheValid) {
            updateVolume();  // Automatically update the volume cache if it's invalid
        }
        return cachedVolume;
    }

    void invalidateVolumeCache() {
        isVolumeCacheValid = false;  // set the cache flag to invalid
    }

};

        // if (T.dimension() != 3) {
        //     std::vector<Point_2> projected_points;
        //     for (const auto& point : points) {
        //         // Project the point to a plane, e.g., the XY-plane
        //         projected_points.emplace_back(point.x(), point.y());
        //     }
        //     Polygon_2 P_2D;
        //     CGAL::convex_hull_2(projected_points.begin(), projected_points.end(), std::back_inserter(P_2D));
        //     // print_2d_convex_hull(P_2D);
        // }
        // else {
        //     try {
        //         CGAL::convex_hull_3_to_face_graph(T, P);
        //         // print_convex_hull(P);
        //     } catch (const std::exception& e) {
        //         std::cerr << "Error computing convex hull: " << e.what() << std::endl;
        //     }
        // }
// try {
//             CGAL::convex_hull_3(points.begin(), points.end(), P);
//         } catch (const std::exception& e) {
//             std::cerr << "Exception: " << e.what() << std::endl;
//             exit(1);
//         }

using ComponentType = std::string;

struct Component {
    std::vector<Triangle> triangles;
    double weight;
    double currentScore;  // Add this line to keep track of the current score
    ComponentType type;  // The type of this component
    int depth = 0;  // Initialize depth to 0 for new components
    Hull convexHull;  // Now it's optional and can be empty
    int id;

    Component() = default;  // Explicitly adding a default constructor

    void addTriangle(const Triangle& triangle) {
        triangles.push_back(triangle);
        updateHull();
        convexHull.invalidateVolumeCache();  // Invalidate the cached volume
    }

    
    // Method to get the last added triangle
    Triangle getLastTriangle() const {
        if (!triangles.empty()) {
            return triangles.back();
        }
        // Handle the case when there are no triangles.
        // You could throw an exception, return a dummy triangle, etc.
        throw std::runtime_error("No triangles in the component");
    }

    void removeTriangle(const Triangle& triangleToRemove) {
        auto it = std::remove_if(triangles.begin(), triangles.end(),
            [&triangleToRemove](const Triangle& triangle) {
                return triangle == triangleToRemove;
            });
        
        if (it != triangles.end()) {
            triangles.erase(it, triangles.end());
            updateHull();
        }
    }

    bool arePointsUnique() const {
        return convexHull.arePointsUnique();
    }

    void updateHull() {
        if (!triangles.empty()) {
            convexHull = Hull(triangles);
            // convexHull.updateVolume();
        }
    }


    bool isValid() const {
        bool contiguous = areTrianglesContiguous(triangles);
        return contiguous;
    }

    void updateWeight() {
        double Zi = static_cast<double>(triangles.size());
        weight = 1.0 / (Zi * Zi);
    }

    std::string toString() const {
        std::ostringstream oss;
        for (const auto& triangle : triangles) {
            oss << triangle.toString() << ", ";
        }
        return oss.str();
    }

    bool overlapsSignificantly(const Component& other);

    // Overload the equality operator
    bool operator==(const Component& other) const {
        // Define your equality logic here. For example:
        return this->triangles == other.triangles; // replace this condition with your specific equality check
    }

    // Overload the inequality operator
    bool operator!=(const Component& other) const {
        return !(*this == other);
    }

    bool isEmpty() const {
        return triangles.empty();
    }

    bool isHullClosed() const {
        return !convexHull.isEmpty() && CGAL::is_closed(convexHull.getSurfaceMesh());
    }


    bool containsTriangle(const Triangle& queryTriangle) const {
        for (const auto& existingTriangle : triangles) {
            if (existingTriangle == queryTriangle) {
                return true;
            }
        }
        return false;
    }

};

struct ComponentHash {
    std::size_t operator()(const Component& component) const {
        std::size_t hashValue = 0;
        for (const auto& triangle : component.triangles) {
            // Assuming Triangle has a hash function
            hashValue ^= std::hash<Triangle>()(triangle) + 0x9e3779b9 + (hashValue << 6) + (hashValue >> 2);
        }
        return hashValue;
    }
};


// Function to identify the component type
ComponentType identifyComponentType(const Component& component) {
    // Check for edge cases
    if (component.triangles.empty()) {
        return "";  // Or some sentinel value to indicate an empty component
    }

    std::unordered_map<std::string, int> frequencyMap;  // Stores the frequency of each triangle type

    // Populate frequencyMap with the frequency of each triangle's type
    for (const auto& triangle : component.triangles) {
        std::string type = triangle.type;  // Use Triangle's type field
        frequencyMap[type]++;
    }

    std::string componentType = "";

    // Construct componentType based on the frequency of each type
    for (const auto& pair : frequencyMap) {
        // Append each type 'pair.first' to the componentType string 'pair.second' times
        componentType += std::string(pair.second, pair.first[0]);  
    }

    // Sort the componentType string to make it standardized
    std::sort(componentType.begin(), componentType.end()); 

    return componentType;
}

bool attemptToMergeComponents(Component& component1, Component& component2, TriangleAdjacency& adjacency) {
    if (&component1 == &component2) {
        return false; // Don't merge the same component with itself
    }

    // Debug: Print triangles in both components before merging
    std::cout << "Triangles in component1 before attempting to merge:\n";
    for (const auto& triangle : component1.triangles) {
        std::cout << triangle.toString() << std::endl;
    }
    std::cout << "Triangles in component2 before attempting to merge:\n";
    for (const auto& triangle : component2.triangles) {
        std::cout << triangle.toString() << std::endl;
    }

    // Create a temporary component with the triangles from both components
    Component tempComponent = component1;
    tempComponent.triangles.insert(tempComponent.triangles.end(), component2.triangles.begin(), component2.triangles.end());

    // Debug: Print triangles in tempComponent after merging
    std::cout << "Triangles in tempComponent after merging:\n";
    for (const auto& triangle : tempComponent.triangles) {
        std::cout << triangle.toString() << std::endl;
    }

    // Check if the merged component is valid
    if (tempComponent.isValid()) {
        // Debug: Print success message
        std::cout << "Merging successful.\n";

        component1 = tempComponent;  // Update component1 to be the merged component
        component1.type = identifyComponentType(component1);  // Update the type of the merged component
        std::cout << "Merged component type: " << component1.type << std::endl;
        return true;
    }

    // Debug: Print failure message
    std::cout << "Merging failed.\n";

    return false;
}


bool Component::overlapsSignificantly(const Component& other) {
    if (this == &other) {
        std::cout << "Comparing the same instance, skipping." << std::endl;
        return false;
    }

    if (*this == other) {
        std::cout << "Comparing identical components based on operator==, skipping." << std::endl;
        return false;
    }

    Hull thisHull = Hull(this->triangles);
    Hull otherHull = Hull(other.triangles);

    double thisVolume = thisHull.volume();
    double otherVolume = otherHull.volume();

    // Handle edge cases where volume is zero
    if (thisVolume == 0.0 || otherVolume == 0.0) {
        std::cout << "Zero volume detected, skipping." << std::endl;
        return false;
    }

    Hull combinedHull(thisHull, otherHull);  // You will need to implement this constructor
    double combinedVolume = combinedHull.volume();

    // Compute the volume that is exclusive to each hull
    double exclusiveThisVolume = combinedVolume - otherVolume;
    double exclusiveOtherVolume = combinedVolume - thisVolume;

    // Check for significant overlap
    double overlapThis = 1.0 - (exclusiveThisVolume / thisVolume);
    double overlapOther = 1.0 - (exclusiveOtherVolume / otherVolume);

    if (overlapThis >= 0.9 || overlapOther >= 0.9) {
        std::cout << "Significant overlap detected: " << overlapThis << ", " << overlapOther << std::endl;
        return true;
    }
    std::cout << "No significant overlap detected: " << overlapThis << ", " << overlapOther << std::endl;
    return false;
}

void mergeOverlappingComponents(std::vector<Component>& components, TriangleAdjacency& adjacency) {
    // exit(1);
    std::set<int> toRemove;  // To keep track of indices to remove

    // Debug: Count total triangles before merging
    int totalTrianglesBefore = 0;
    for (const auto& component : components) {
        totalTrianglesBefore += component.triangles.size();
    }
    std::cout << "Total triangles before merging: " << totalTrianglesBefore << std::endl;

    // Store volumes to avoid redundant calculations
    std::vector<double> volumes(components.size(), 0.0);
    for (int i = 0; i < components.size(); ++i) {
        Hull hull_i = Hull(components[i].triangles);
        volumes[i] = hull_i.volume();
    }

    for (int i = 0; i < components.size(); ++i) {
        if (toRemove.count(i) > 0) {
            continue;  // Skip already removed components
        }

        double volume_i = volumes[i];

        for (int j = i + 1; j < components.size(); ++j) {
            if (toRemove.count(j) > 0) {
                continue;  // Skip already removed components
            }

            double volume_j = volumes[j];

            // Check for type and significant overlap
            if (components[i].type == components[j].type && components[i].overlapsSignificantly(components[j])) {
                if (attemptToMergeComponents(components[i], components[j], adjacency)) {
                    toRemove.insert(j);  // Mark for removal
                }
            }
        }
    }

    // Remove marked components
    for (auto it = toRemove.rbegin(); it != toRemove.rend(); ++it) {
        components.erase(components.begin() + *it);
    }

    // Debug: Count total triangles after merging
    int totalTrianglesAfter = 0;
    for (const auto& component : components) {
        totalTrianglesAfter += component.triangles.size();
    }
    std::cout << "Total triangles after merging: " << totalTrianglesAfter << std::endl;

    if (totalTrianglesBefore != totalTrianglesAfter) {
        std::cout << "Warning: Some triangles are missing!" << std::endl;
    }
}



double ATaTb(const Cluster& Ta, const Cluster& Tb, 
             const std::unordered_map<Triangle, std::vector<Triangle>>& potentialNeighborsCache, 
             const TriangleAdjacency& triangleAdjacency,
             const std::unordered_set<Triangle>& tbTriangles,
             double tau_N) {  // <-- Added tau_N for early exit
    int count = 0;
    int minSize = std::min(Ta.triangles.size(), Tb.triangles.size());
    double earlyExitThreshold = tau_N * minSize;  // Calculate the threshold for early exit

    for (const auto& ta_i : Ta.triangles) {
        auto potentialNeighbors = potentialNeighborsCache.at(ta_i);
        for (const auto& potentialNeighbor : potentialNeighbors) {
            if (tbTriangles.count(potentialNeighbor) && ta_i.isAdjacent(potentialNeighbor, triangleAdjacency)) {
                count++;
                // if (count >= earlyExitThreshold) {
                //     return 1.0;  // or any value >= tau_N to indicate high similarity
                // }
            }
        }
    }

    double result = static_cast<double>(count) / minSize;
    return result;
}



void saveToOBJ(const std::vector<Cluster>& clusters, const std::string& filename) {
    std::ostringstream objStream;
    std::ostringstream mtlStream;

    // Write MTL file
    std::string mtlFilename = filename + ".mtl";
    objStream << "mtllib " << mtlFilename << "\n";

    int materialIndex = 0;
    size_t totalVerticesProcessed = 0;  // Track total vertices processed so far

    for (const auto& cluster : clusters) {
        if (cluster.triangles.empty()) continue;

        // Assume the first triangle of each cluster holds the color for the entire cluster
        Vector3D color = cluster.triangles[0].color;
        mtlStream << "newmtl material" << materialIndex << "\n";
        mtlStream << "Kd " << color.x << " " << color.y << " " << color.z << "\n\n";
        
        for (const auto& triangle : cluster.triangles) {
            objStream << "v " << triangle.vec1.x << " " << triangle.vec1.y << " " << triangle.vec1.z << "\n";
            objStream << "v " << triangle.vec2.x << " " << triangle.vec2.y << " " << triangle.vec2.z << "\n";
            objStream << "v " << triangle.vec3.x << " " << triangle.vec3.y << " " << triangle.vec3.z << "\n";
        }

        objStream << "usemtl material" << materialIndex << "\n";
        for (size_t j = 0; j < cluster.triangles.size(); ++j) {
            size_t baseIdx = 1 + j * 3 + totalVerticesProcessed;  // Adjusted index calculation
            objStream << "f " << baseIdx << " " << (baseIdx + 1) << " " << (baseIdx + 2) << "\n";
        }

        totalVerticesProcessed += cluster.triangles.size() * 3;  // Update the total vertex count
        ++materialIndex;
    }

    // Save to files
    std::ofstream objFile(filename + ".obj");
    objFile << objStream.str();
    objFile.close();

    std::ofstream mtlFile(mtlFilename);
    mtlFile << mtlStream.str();
    mtlFile.close();
}

void saveComponentsToOBJ(const std::vector<Component>& components, const std::string& filename) {
    std::ostringstream objStream;
    std::ostringstream mtlStream;

    // Write MTL file
    std::string mtlFilename = filename + ".mtl";
    objStream << "mtllib " << mtlFilename << "\n";

    int materialIndex = 0;
    size_t totalVerticesProcessed = 0;  // Track total vertices processed so far

    for (const auto& component : components) {
        if (component.triangles.empty()) continue;

        // Assume the first triangle of each component holds the color for the entire component
        Vector3D color = component.triangles[0].color;
        mtlStream << "newmtl material" << materialIndex << "\n";
        mtlStream << "Kd " << color.x << " " << color.y << " " << color.z << "\n\n";
        
        for (const auto& triangle : component.triangles) {
            objStream << "v " << triangle.vec1.x << " " << triangle.vec1.y << " " << triangle.vec1.z << "\n";
            objStream << "v " << triangle.vec2.x << " " << triangle.vec2.y << " " << triangle.vec2.z << "\n";
            objStream << "v " << triangle.vec3.x << " " << triangle.vec3.y << " " << triangle.vec3.z << "\n";
        }

        objStream << "usemtl material" << materialIndex << "\n";
        for (size_t j = 0; j < component.triangles.size(); ++j) {
            size_t baseIdx = 1 + j * 3 + totalVerticesProcessed;  // Adjusted index calculation
            objStream << "f " << baseIdx << " " << (baseIdx + 1) << " " << (baseIdx + 2) << "\n";
        }

        totalVerticesProcessed += component.triangles.size() * 3;  // Update the total vertex count
        ++materialIndex;
    }

    // Save to files
    std::ofstream objFile(filename + ".obj");
    objFile << objStream.str();
    objFile.close();

    std::ofstream mtlFile(mtlFilename);
    mtlFile << mtlStream.str();
    mtlFile.close();
}

double minMaxScale(double x, double min_val, double max_val) {
    if (max_val - min_val == 0) return 0;
    return (x - min_val) / (max_val - min_val);
}


double Stitj(const Triangle& ti, const Triangle& tj, const std::vector<double> weights) {
    double wA = weights[0], wL = weights[1], wS = weights[2], wN = weights[3];
    
    double areaDenominator = std::max(ti.area, tj.area);
    double maxLengthDenominator = ti.e_max + tj.e_max;
    double minLengthDenominator = ti.e_min + tj.e_min;

    double areaDifference = (areaDenominator != 0.0) ? wA * std::abs(ti.area - tj.area) / areaDenominator : 0.0;
    double maxLengthDifference = (maxLengthDenominator != 0.0) ? wL * std::abs(ti.e_max - tj.e_max) / maxLengthDenominator : 0.0;
    double minLengthDifference = (minLengthDenominator != 0.0) ? wS * std::abs(ti.e_min - tj.e_min) / minLengthDenominator : 0.0;
    // Normalize the normals
    Vector3D normal1 = ti.normal;
    Vector3D normal2 = tj.normal;
    normal1.normalize();
    normal2.normalize();
    double dotProduct = normal1.dot(normal2);
    
    double normalDifference = wN * (1 - dotProduct) / 2.0;

    double totalSum = areaDifference + maxLengthDifference + minLengthDifference + normalDifference;

    // std::cout << "Normal vector of ti: (" << ti.normal.x << ", " << ti.normal.y << ", " << ti.normal.z << ")\n";
    // std::cout << "Normal vector of tj: (" << tj.normal.x << ", " << tj.normal.y << ", " << tj.normal.z << ")\n";
    // std::cout << "Dot product: " << ti.normal.dot(tj.normal) << "\n";

    return totalSum;
}

void writeVectorToFile(const std::vector<double>& vec, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Could not open file: " << filename << std::endl;
        return;
    }

    for (const auto& elem : vec) {
        outFile << elem << "\n";
    }

    outFile.close();
}

std::vector<Cluster> initialClusteringByShapeSimilarity(const std::vector<Triangle>& triangles, std::vector<double> weights, double tau_S) {
    // Sort triangles in order of increasing area
    std::vector<Triangle> sortedTriangles = triangles;
    std::sort(sortedTriangles.begin(), sortedTriangles.end(), [](const Triangle& a, const Triangle& b) {
        return a.area < b.area;
    });

    std::vector<Cluster> clusters;

    for (size_t i = 0; i < sortedTriangles.size(); ++i) {
        bool assigned = false;
        for (Cluster& cluster : clusters) {
            for (const Triangle& existingTriangle : cluster.triangles) {  // Compare with every triangle in the cluster
                double similarity = Stitj(sortedTriangles[i], existingTriangle, weights);
                if (similarity <= tau_S) {
                    cluster.triangles.emplace_back(sortedTriangles[i]);
                    assigned = true;
                    break;  // No need to check other triangles in this cluster
                }
            }
            if (assigned) {
                break;  // No need to check other clusters
            }
        }
        if (!assigned) {
            // Create a new cluster for this triangle
            clusters.emplace_back(Cluster{{sortedTriangles[i]}});
        }
    }
    // writeVectorToFile(all_stitj, "all_stitj_non_normalize.txt"); 
    return clusters;
}

struct BacktrackingState {
    int nextComponentToBacktrack = 0;
    int individualRotationsDone = 0;
    bool multiComponentRotation = false;
    std::unordered_set<std::string> triedCombinations;
    std::map<int, std::deque<Triangle>> trianglesTriedInThisState;
    std::set<int> componentsAlreadyTried;
};

struct AlgorithmState {
    std::vector<Component> tempComponents;
    std::unordered_set<Triangle> remainingTriangles;
    std::unordered_set<Triangle> usedTriangles;
    std::vector<int> repetitionCounts;
    std::vector<double> componentWeights;
    std::vector<double> globalVolumeChangeHistory;

    std::vector<std::vector<std::pair<double, Triangle>>> adjustedCandidateTriangles;
    std::vector<std::stack<std::vector<Triangle>>> backtrackCandidateStacks;
    std::vector<bool> componentsStalled;
    std::stack<int> componentIndexStack;
    std::vector<std::stack<Triangle>> issueTrianglesStacks;
    std::unordered_set<std::string> triedCombinations;
    std::map<int, std::deque<Triangle>> trianglesThatLedToThisState;
    std::vector<std::vector<std::pair<Triangle, double>>> candidateTriangles;
    bool anyComponentGrown;
    bool allComponentsStalled;

    int nextComponentToBacktrack = 0; // add this new field
};

void printSet(const std::unordered_set<Triangle>& mySet) {
    std::cout << "{";
    for (auto it = mySet.begin(); it != mySet.end(); ++it) {
        if (it != mySet.begin()) {
            std::cout << ", ";
        }
        std::cout << it->toString();
    }
    std::cout << "}" << std::endl;
}

void printCandidateTrianglesForComponents(const std::vector<std::vector<std::pair<Triangle, double>>>& candidateTrianglesForComponents) {
    for (size_t i = 0; i < candidateTrianglesForComponents.size(); ++i) {
        std::cout << "Component " << i<< ":\n";
        for (const auto& pair : candidateTrianglesForComponents[i]) {
            Triangle triangle = pair.first;
            double score = pair.second;
            std::cout << "Triangle: " << triangle.toString() << ", Score: " << score << "\n";
        }
        std::cout << std::endl;
    }
}

// New function to compute volume changes for all triangles and a given component
std::vector<std::pair<double, Triangle>> computeVolumeChanges(
    const Component& component, 
    const std::unordered_set<Triangle>& remainingTriangles
) {
    std::vector<std::pair<double, Triangle>> volumeChanges;
    double originalVolume = component.convexHull.volume();

    Component tempComponent = component;  // Local, non-const copy
    for (const Triangle& triangle : remainingTriangles) {
        tempComponent.addTriangle(triangle);  // Work with local copy
        
        if (!tempComponent.isValid()) {
            tempComponent.removeTriangle(triangle);
            continue;
        }
        
        double newVolume = tempComponent.convexHull.volume();
        double volumeChange = std::abs(newVolume - originalVolume);
        
        if (volumeChange != 0) {
            volumeChanges.emplace_back(volumeChange, triangle);
        }
        
        tempComponent.removeTriangle(triangle);
    }
    
    return volumeChanges;
}



std::vector<std::vector<std::pair<Triangle, double>>> getCandidateTrianglesForEachComponent(
    const AlgorithmState& currentState,
    size_t topN = 2
) {
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<std::pair<Triangle, double>>> candidateTrianglesForComponents;

    std::unordered_set<Triangle> usedTriangles = currentState.usedTriangles; // Assuming usedTriangles is hashable

    auto preloop = std::chrono::high_resolution_clock::now();
    for (const Component& component : currentState.tempComponents) {
        std::priority_queue<std::pair<double, Triangle>> candidatePriorityQueue;

        auto volumeStart = std::chrono::high_resolution_clock::now();
        auto volumeChanges = computeVolumeChanges(component, currentState.remainingTriangles);  // Pass the const reference directly
        auto volumeEnd = std::chrono::high_resolution_clock::now();
        std::cout << "Volume loop time: " << std::chrono::duration_cast<std::chrono::milliseconds>(volumeEnd - volumeStart).count() << " ms" << std::endl;
        std::sort(volumeChanges.begin(), volumeChanges.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });
        std::vector<std::pair<Triangle, double>> topCandidates;
        for (size_t i = 0; i < topN && i < volumeChanges.size(); ++i) {
            topCandidates.push_back({volumeChanges[i].second, volumeChanges[i].first});
        }
        candidateTrianglesForComponents.push_back(topCandidates);
        
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Total time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

    return candidateTrianglesForComponents;
}


double calculateConvexHullChange(Component& component, const Triangle& candidateTriangle) {
    double originalVolume = component.convexHull.volume();

    // Create a new Hull object by first copying the existing component
    Component tempComponent = component;

    // Then add the new triangle to this temporary component
    tempComponent.addTriangle(candidateTriangle);  // Assuming you have a way to add a triangle to a Component

    // Compute the hull for this new component
    // (automatically done within addTriangle())
    
    double newVolume = tempComponent.convexHull.volume();
    return std::abs(newVolume - originalVolume);  // Returns absolute change in volume
}


// Define a struct to hold a set of best triangles along with their global change value
struct BestTriangleSet {
    double globalChange;
    std::vector<Triangle> triangles;
    
    bool operator<(const BestTriangleSet& other) const {
        return globalChange < other.globalChange;
    }
};

void findBestTrianglesRecursively(
    const std::vector<std::vector<std::pair<double, Triangle>>>& candidatesForComponents,
    std::vector<Triangle>& currentCombination,
    double currentGlobalChange,
    double globalRepetitionFactor,
    std::priority_queue<BestTriangleSet>& bestTriangleSets,
    size_t componentIndex = 0) {
    
    if (componentIndex == candidatesForComponents.size()) {
        // Base case: We've selected one triangle for each component
        BestTriangleSet newSet;
        newSet.globalChange = currentGlobalChange * globalRepetitionFactor;
        newSet.triangles = currentCombination;
        bestTriangleSets.push(newSet);
        return;
    }
    
    for (const auto& candidate : candidatesForComponents[componentIndex]) {
        currentCombination[componentIndex] = candidate.second;
        findBestTrianglesRecursively(
            candidatesForComponents,
            currentCombination,
            currentGlobalChange + candidate.first,
            globalRepetitionFactor,
            bestTriangleSets,
            componentIndex + 1
        );
    }
}
// void findBestTrianglesRecursively(
//     const std::vector<std::vector<std::pair<double, Triangle>>>& candidatesForComponents,
//     std::vector<Triangle>& currentCombination,
//     double currentGlobalChange,
//     double globalRepetitionFactor,
//     std::unordered_map<std::string, int>& componentFrequencies,  // Key: component type identifier, Value: frequency of repetition
//     std::priority_queue<BestTriangleSet>& bestTriangleSets,
//     size_t componentIndex = 0) {
    
//     if (componentIndex == candidatesForComponents.size()) {
//         // Base case: We've selected one triangle for each component
//         double finalGlobalChange = currentGlobalChange * globalRepetitionFactor;
        
//         for (size_t i = 0; i < currentCombination.size(); ++i) {
//             std::string componentIdentifier = /* some way to identify the component */;
            
//             int frequency = componentFrequencies[componentIdentifier];  // Fetch the frequency of the component
//             double weight = 1.0 / (std::pow(frequency, 2));  // Compute the weight according to Eqn. (4)
            
//             finalGlobalChange *= weight;  // Update the global change using the weight
//         }
        
//         BestTriangleSet newSet;
//         newSet.globalChange = finalGlobalChange;
//         newSet.triangles = currentCombination;
//         bestTriangleSets.push(newSet);
        
//         // Update componentFrequencies based on this new successful combination
//         // ...
        
//         return;
//     }
    
//     for (const auto& candidate : candidatesForComponents[componentIndex]) {
//         currentCombination[componentIndex] = candidate.second;
//         findBestTrianglesRecursively(
//             candidatesForComponents,
//             currentCombination,
//             currentGlobalChange + candidate.first,
//             globalRepetitionFactor,
//             componentFrequencies,
//             bestTriangleSets,
//             componentIndex + 1
//         );
//     }
// }


static bool isValidAfterAddingTriangle(const Component& component, const Triangle& newTriangle) {
    // Check if the triangle is already in the component
    if (std::find(component.triangles.begin(), component.triangles.end(), newTriangle) != component.triangles.end()) {
        return false;
    }

    // Create a temporary list of triangles that includes the new triangle
    std::vector<Triangle> tempTriangles = component.triangles;
    tempTriangles.push_back(newTriangle);
    return areTrianglesContiguous(tempTriangles);  // Assuming you can call this function directly
}


double calculateGlobalRepetitionFactor(
    const std::vector<double>& globalVolumeChangeHistory,
    const std::vector<double>& minConvexHullChanges) {
    
    // If there's not enough data, return a neutral factor
    if (globalVolumeChangeHistory.empty() || minConvexHullChanges.empty()) {
        return 1.0;
    }

    // Calculate the mean and variance of the historical volume changes
    double mean = std::accumulate(globalVolumeChangeHistory.begin(), globalVolumeChangeHistory.end(), 0.0) / globalVolumeChangeHistory.size();
    double sq_sum = std::inner_product(globalVolumeChangeHistory.begin(), globalVolumeChangeHistory.end(), globalVolumeChangeHistory.begin(), 0.0);
    double variance = (sq_sum / globalVolumeChangeHistory.size()) - (mean * mean);

    // If variance is zero, the system is perfectly repetitive
    if (variance == 0) {
        return 1.0;
    }

    // Calculate the mean of the current minimum changes
    double currentMean = std::accumulate(minConvexHullChanges.begin(), minConvexHullChanges.end(), 0.0) / minConvexHullChanges.size();
  
    // Calculate a repetition factor based on how close the current mean is to the historical mean, and how low the variance is
    double distanceToMean = std::abs(currentMean - mean);
    double repetitionFactor = std::exp(-distanceToMean / variance);

    return repetitionFactor;
}

std::vector<std::vector<Triangle>> synchronizeGrowth(
    AlgorithmState& currentState,
    const std::vector<double>& minConvexHullChanges,
    bool considerRepetition) {

    // Directly use currentState.adjustedCandidateTriangles without copying
    auto& candidatesForComponents = currentState.adjustedCandidateTriangles;

    double globalRepetitionFactor = 1.0;
    if (considerRepetition) {
        globalRepetitionFactor = calculateGlobalRepetitionFactor(currentState.globalVolumeChangeHistory, minConvexHullChanges);
    }

    std::priority_queue<BestTriangleSet> bestTriangleSets;
    std::vector<Triangle> currentCombination(currentState.tempComponents.size());

    findBestTrianglesRecursively(
        candidatesForComponents,
        currentCombination,
        0.0,
        globalRepetitionFactor,
        bestTriangleSets
    );

    std::vector<std::unordered_set<Triangle>> uniqueTrianglesForComponents(currentState.tempComponents.size());
    std::vector<std::vector<Triangle>> bestTrianglesForComponents(currentState.tempComponents.size());

    while (!bestTriangleSets.empty()) {
        const BestTriangleSet& bestSet = bestTriangleSets.top();
        for (size_t i = 0; i < currentState.tempComponents.size(); ++i) {
            if (uniqueTrianglesForComponents[i].insert(bestSet.triangles[i]).second) {
                bestTrianglesForComponents[i].push_back(bestSet.triangles[i]);
            }
        }
        bestTriangleSets.pop();
    }
    return bestTrianglesForComponents;
}


std::vector<Triangle> findLargestDisjointSimilarTriangles(const std::vector<Triangle>& triangles, std::vector<Triangle>& seeds, 
                                                          std::vector<double> weights, double tau_S) {

    if (triangles.empty()) {
        std::cerr << "Error: No triangles to find seeds from.\n";
        return seeds;
    }

    // Sort triangles by area in descending order
    std::vector<Triangle> sortedTriangles = triangles;
    std::sort(sortedTriangles.begin(), sortedTriangles.end(), [](const Triangle& a, const Triangle& b) {
        return a.area > b.area;
    });

    // std::cout << "weight is: " << weights[0] << " " << weights[1] << " " << weights[2] << " " << weights[3] << std::endl;
    // exit(1);

    for (const Triangle& triangle : sortedTriangles) {
        bool isSimilar = true;
        for (const Triangle& seed : seeds) {
            double similarity = Stitj(triangle, seed, weights);
            if (similarity > tau_S) { 
                isSimilar = false;
                break;
            }
        }
        if (isSimilar) {
            seeds.push_back(triangle);
        }
    }

    // If D is empty, add the largest triangle as the only seed
    if (seeds.empty()) {
        seeds.push_back(sortedTriangles.front());
    }

    return seeds;
}



std::vector<Triangle> findSeeds(const std::unordered_set<Triangle>& remainingTriangles,
                                std::vector<Triangle>& seeds,
                                std::vector<double> weights,
                                double tau_S) {
    std::cout << "Inside findSeeds, remainingTriangles size: " << remainingTriangles.size() << std::endl;
    
    // Debug: Print remainingTriangles
    // std::cout << "Remaining Triangles:\n";
    // for (const auto& tri : remainingTriangles) {
    //     std::cout << tri.toString() << "\n";  // Assumes Triangle has an overloaded << operator
    // }

    auto result = findLargestDisjointSimilarTriangles(
        std::vector<Triangle>(remainingTriangles.begin(), remainingTriangles.end()),
        seeds,
        weights,
        tau_S
    );
    
    // Debug: Print seeds
    std::cout << "Seeds:\n";
    for (const auto& seed : seeds) {
        std::cout << seed.toString() << "\n";  // Assumes Triangle has an overloaded << operator
    }

    std::cout << "Number of seeds found: " << result.size() << std::endl;
    return result;
}



// std::vector<Component> initializeTempComponents(const std::vector<Triangle>& seeds,
//                                                 std::unordered_set<Triangle>& remainingTriangles) {
//     std::vector<Component> tempComponents;

//     if (seeds.size() == 1) {
//         std::cout << "Warning: Only one seed found. This may affect the performance or correctness of the algorithm." << std::endl;
//         // exit(1);
//         // Optionally, you can also choose to return an empty vector here or continue execution.
//     }

//     for (const auto& seed : seeds) {
//         auto eraseCount = remainingTriangles.erase(seed);
//         if (eraseCount != 0) {
//             Component newComponent;  // Should work now with the default constructor
//             newComponent.addTriangle(seed);  // This will also update the Hull
//             tempComponents.push_back(newComponent);

//             // Debug: Successfully added triangle to component
//             std::cout << "Debug: Added triangle to new component: " << seed.toString() << std::endl;
//         } else {
//             // Debug: Triangle was not in remainingTriangles
//             std::cout << "Debug: Seed was not in remainingTriangles: " << seed.toString() << std::endl;
//         }
//     }

//     // Debug: Number of remaining triangles after initialization
//     std::cout << "Debug: Number of remainingTriangles after initialization: " << remainingTriangles.size() << std::endl;

//     return tempComponents;
// }

std::vector<Component> initializeTempComponents(const std::vector<Triangle>& seeds,
                                                std::unordered_set<Triangle>& remainingTriangles) {
    std::cout << "Inside initializeTempComponents, remainingTriangles size: " << remainingTriangles.size() << std::endl;
    std::vector<Component> tempComponents;
    for (const auto& seed : seeds) {
        auto eraseCount = remainingTriangles.erase(seed);
        std::cout << "erased seed from remaining triangle " << seed.toString() << std::endl;
        if (eraseCount != 0) {
            Component newComponent;
            newComponent.addTriangle(seed);
            tempComponents.push_back(newComponent);
        }
    }
    std::cout << "after initializeTempComponents, remainingTriangles size: " << remainingTriangles.size() << std::endl;
    return tempComponents;
}



double calculateWeight(int repetitionCount) {
    return 1.0 / (std::pow(repetitionCount, 2));
}


bool isValidWeight(double newWi, double oldWi) {

    // We only check for less than or equal to because lower weight is better
    // and the weight may not change (stay equal) in some iterations
    return (newWi <= oldWi);
}


std::string serializePathVector(const std::vector<Component>& components) {
    std::stringstream ss;
    for (size_t i = 0; i < components.size(); ++i) {
        ss << components[i].toString();
        if (i < components.size() - 1) {
            ss << "|";
        }
    }
    return ss.str();
}

void inspectStack(std::stack<std::vector<Component>> &originalStack) {
    std::stack<std::vector<Component>> tempStack;

    while (!originalStack.empty()) {
        auto &topElement = originalStack.top();
        
        std::cout << "Vector size: " << topElement.size() << std::endl;

        for (size_t i = 0; i < topElement.size(); ++i) {
            const Component& comp = topElement[i];
            // std::cout << "Component ID: " << comp.id << std::endl;  // Access the ID here
            std::cout << "Component " << comp.id << " contains the following triangles:" << std::endl;

            for (const auto& triangle : comp.triangles) {
                std::cout << "Triangle: " << triangle.toString() << std::endl;
            }
        }

        tempStack.push(topElement);
        originalStack.pop();
    }

    // Restore original stack
    while (!tempStack.empty()) {
        originalStack.push(tempStack.top());
        tempStack.pop();
    }
}


void printStack(std::stack<Triangle> s) {
    while (!s.empty()) {
        std::cout << s.top().toString() << std::endl;
        s.pop();  // Note: This will empty the stack
    }
    std::cout << "----" << std::endl;
}

void printIndexStack(std::stack<int> &originalStack) {
    std::stack<int> tempStack;

    // Print elements and push them into a temporary stack
    while (!originalStack.empty()) {
        int topElement = originalStack.top();
        std::cout << "Component Index: " << topElement << std::endl;
        tempStack.push(topElement);
        originalStack.pop();
    }

    // Restore the original stack from the temporary stack
    while (!tempStack.empty()) {
        originalStack.push(tempStack.top());
        tempStack.pop();
    }
}

void printBestTrianglesForComponents(const std::vector<std::vector<Triangle>>& bestTrianglesForComponents) {
    for (size_t i = 0; i < bestTrianglesForComponents.size(); ++i) {
        std::cout << "Component " << i << ":\n";
        for (const Triangle& triangle : bestTrianglesForComponents[i]) {
            std::cout << "Triangle: " << triangle.toString() << "\n";
        }
        std::cout << std::endl;
    }
}

std::string generateKeyForCombination(const std::map<int, std::deque<Triangle>> &trianglesMap) {
    std::string key;
    for (const auto& pair : trianglesMap) {
        int componentIndex = pair.first;
        const std::deque<Triangle>& triangles = pair.second;
        for (const auto& triangle : triangles) {
            key += triangle.toString() + ";";
        }
    }
    return key;
}

int countOccurrences(const std::string& str, const std::string& sub) {
    int count = 0;
    size_t pos = 0;
    while ((pos = str.find(sub, pos)) != std::string::npos) {
        ++count;
        pos += sub.size();
    }
    return count;
}


void printGroupedByComponent(const std::vector<std::vector<Triangle>>& groupedByComponent) {
    for (size_t i = 0; i < groupedByComponent.size(); ++i) {
        std::cout << "Component " << i << ": ";
        
        for (size_t j = 0; j < groupedByComponent[i].size(); ++j) {
            // Assume Triangle has a method toString() to return its string representation
            std::cout << groupedByComponent[i][j].toString();
            
            // Print comma if not the last element
            if (j < groupedByComponent[i].size() - 1) {
                std::cout << ", ";
            }
        }
        
        std::cout << std::endl;
    }
}

// bool handleRotatingComponent(AlgorithmState& currentState,
//                              int componentToBacktrack,
//                              std::vector<int>& triangleIndices,
//                              std::vector<std::pair<int, Triangle>>& currentCombinationTriangles,
//                              BacktrackingState& backtrackingState) {

//     // Keep track of components that can be rotated
//     std::vector<int> rotatableComponents;
    
//     for (size_t i = 0; i < currentState.backtrackCandidateStacks.size(); ++i) {
//         if (currentState.backtrackCandidateStacks[i].size() >= 2) {
//             rotatableComponents.push_back(i);
//         }
//     }

//     // If no component can be rotated, it's a failure
//     if (rotatableComponents.empty()) {
//         std::cout << "No component can be rotated" << std::endl;
//         return false;
//     }
    
//     for (size_t i = 0; i < currentState.backtrackCandidateStacks.size(); ++i) {
//         std::vector<Triangle> componentTriangles = currentState.backtrackCandidateStacks[i].top();
//         Triangle selectedTriangle;

//         // Check if this component can be rotated
//         if (std::find(rotatableComponents.begin(), rotatableComponents.end(), i) != rotatableComponents.end()) {
//             if (backtrackingState.multiComponentRotation || i == componentToBacktrack) {
//                 selectedTriangle = componentTriangles[1];
//             } else {
//                 selectedTriangle = componentTriangles[0];
//             }
//         } else {
//             if (componentTriangles.empty()) {
//                 continue; // Skip this component, as it can't contribute a triangle
//             }
//             selectedTriangle = componentTriangles[0]; // Default option
//         }

//         currentCombinationTriangles.emplace_back(i, selectedTriangle);
//     }

//     return true;
// }
bool handleRotatingComponent(AlgorithmState& currentState,
                             int componentToBacktrack,
                             std::vector<int>& triangleIndices,
                             std::vector<std::pair<int, Triangle>>& currentCombinationTriangles,
                             BacktrackingState& backtrackingState) {
    
    std::cout << "Handling component " << componentToBacktrack << " for rotation" << std::endl;  // Debugging line

    for (size_t i = 0; i < currentState.backtrackCandidateStacks.size(); ++i) {
        if (currentState.backtrackCandidateStacks[i].empty()) {
            continue;
        }

        std::vector<Triangle> componentTriangles = currentState.backtrackCandidateStacks[i].top();
        Triangle selectedTriangle;

        if (backtrackingState.multiComponentRotation) {
            int indexToRotate = 1 % componentTriangles.size();  // Always rotate to the second triangle if it exists.
            selectedTriangle = componentTriangles[indexToRotate];
        } else {
            if (i == componentToBacktrack && componentTriangles.size() >= 2) {
                selectedTriangle = componentTriangles[1];  // Single-component rotation for components with multiple triangles
            } else {
                selectedTriangle = componentTriangles[0];  // Default option (also for components with only 1 triangle)
            }
        }

        currentCombinationTriangles.emplace_back(i, selectedTriangle);
    }

    return true;
}



bool recursiveBackTracking(std::stack<AlgorithmState>& globalBacktrackStack,
                            AlgorithmState& currentState,
                            BacktrackingState &backtrackingState) {
    if (globalBacktrackStack.empty()) {
        std::cout << "Global backtrack stack is empty" << std::endl;
        exit(1);
        return false;  // No solution exists
    }
    
    std::vector<int> triangleIndices(currentState.tempComponents.size(), 0);

    while (true) {
        // std::cout << "Current rotating component: " << backtrackingState.nextComponentToBacktrack << std::endl;  // Debug line 1
        std::cout << "STARTING LOOP " << std::endl;
        std::vector<std::pair<int, Triangle>> currentCombinationTriangles;
        if (!handleRotatingComponent(currentState, backtrackingState.nextComponentToBacktrack, triangleIndices, currentCombinationTriangles, backtrackingState)) {
            // Record the failed state if you need to.
            std::string failedKey = generateKeyForCombination(backtrackingState.trianglesTriedInThisState);
            backtrackingState.triedCombinations.insert(failedKey);
            // Clear currentCombinationTriangles to keep it in a clean state.
            currentCombinationTriangles.clear();
            std::cout << "Failed to rotate component " << backtrackingState.nextComponentToBacktrack << std::endl;
            exit(1);
            return false;
        }
        std::vector<std::vector<Triangle>> groupedByComponent(currentState.tempComponents.size());

        for (const auto& pair : currentCombinationTriangles) {
            int componentIndex = pair.first;
            Triangle triangle = pair.second;
            groupedByComponent[componentIndex].push_back(triangle);
            backtrackingState.trianglesTriedInThisState[componentIndex].push_back(triangle);
        }
        std::string key = generateKeyForCombination(backtrackingState.trianglesTriedInThisState); // Adjusted to use BacktrackingState
        if (backtrackingState.triedCombinations.find(key) != backtrackingState.triedCombinations.end()) {
            std::cout << "Already visited this combination: " << key << std::endl;
            // exit(1);
            return false;
        } else {
            std::cout << "This is a new combination: " << key << std::endl; // Debugging line
        }
        backtrackingState.triedCombinations.insert(key);
        // Clearing trianglesTriedInThisState
        for(auto& pair : backtrackingState.trianglesTriedInThisState) {
            pair.second.clear();
        }
        // printGroupedByComponent(groupedByComponent);
        bool successfullyAddedTriangles = false;
        for (size_t i = 0; i < currentState.tempComponents.size(); ++i) {
            if (groupedByComponent[i].empty()) {
                std::cout << "EMPTYYY" << std::endl;
                continue;
            }
            // Initialize vectors for hypothetical updates
            std::vector<int> hypotheticalNewZi(currentState.tempComponents.size());
            std::vector<double> hypotheticalNewWi(currentState.tempComponents.size());

            // Calculate hypothetical new values
            for (size_t j = 0; j < currentState.tempComponents.size(); ++j) {
                int oldZi = currentState.repetitionCounts[j];
                hypotheticalNewZi[j] = oldZi + 1;
                double oldWi = currentState.componentWeights[j];
                hypotheticalNewWi[j] = calculateWeight(hypotheticalNewZi[j]);
            }

            // Validate triangles and weights
            bool allUpdatesAreValid = true;
            bool allBestTrianglesAreValid = true;
            for (size_t j = 0; j < currentState.tempComponents.size(); ++j) {
                Triangle bestTriangleForValidation = groupedByComponent[j][0]; // Assumes groupedByComponent[i] is not empty

                bool validTriangle = isValidAfterAddingTriangle(currentState.tempComponents[j], bestTriangleForValidation);
                bool validWeight = isValidWeight(hypotheticalNewWi[j], currentState.componentWeights[j]);

                if (!validTriangle || !validWeight) {
                    std::cout << "Triangle " << bestTriangleForValidation.toString() << " is not valid" << std::endl;

                    backtrackingState.nextComponentToBacktrack = (backtrackingState.nextComponentToBacktrack + 1) % currentState.tempComponents.size();
                    backtrackingState.individualRotationsDone++;

                    if (backtrackingState.individualRotationsDone >= currentState.tempComponents.size()) {
                        std::cout << "All rotations done" << std::endl;
                        backtrackingState.multiComponentRotation = true;
                    }

                    return recursiveBackTracking(globalBacktrackStack, currentState, backtrackingState);
                }

                if (currentState.usedTriangles.find(bestTriangleForValidation) != currentState.usedTriangles.end()) {
                    std::cout << "Triangle " << bestTriangleForValidation.toString() << " is already used" << std::endl;
                    // You may want to handle this case, perhaps setting a flag or initiating another round of backtracking.
                }
            }


            if (allUpdatesAreValid && allBestTrianglesAreValid) {
                for (size_t j = 0; j < currentState.tempComponents.size(); ++j) {
                    // Assuming validation passed, add the triangle to the component
                    successfullyAddedTriangles = true;
                    Triangle bestTriangleForValidation = groupedByComponent[j][0];
                    currentState.tempComponents[j].addTriangle(bestTriangleForValidation);
                    currentState.trianglesThatLedToThisState[j].push_back(bestTriangleForValidation);

                    // Update other state variables
                    currentState.repetitionCounts[j] = hypotheticalNewZi[j];
                    currentState.componentWeights[j] = hypotheticalNewWi[j];
                    currentState.usedTriangles.insert(bestTriangleForValidation);
                    currentState.remainingTriangles.erase(bestTriangleForValidation);
                    std::cout << "all updates are valid here" << std::endl;
                    std::cout << "best triangle is " << bestTriangleForValidation.toString() << std::endl;
                }
                globalBacktrackStack.push(currentState);
                currentState.candidateTriangles = getCandidateTrianglesForEachComponent(currentState);
                printCandidateTrianglesForComponents(currentState.candidateTriangles);
                bool allComponentsHaveCandidates = true;
                for (const auto& candidates : currentState.candidateTriangles) {
                    if (candidates.empty()) {
                        allComponentsHaveCandidates = false;
                        globalBacktrackStack.pop();
                        if (globalBacktrackStack.empty()) {
                            std::cout << "The global backtrack stack is empty. Cannot proceed." << std::endl;
                            exit(1);
                            return false;
                        }
                        currentState = globalBacktrackStack.top();
                        return false;
                    }
                }
                if (allComponentsHaveCandidates && successfullyAddedTriangles) {
                    return true;
                }
            }
            //Need to fix this and recursively go back maybe? or add this key to triedCombinations 
            else {
                std::cout << "Updates are not valid. Reverting to previous state." << std::endl;
                exit(1);
            }
        }
    }
}

// Custom comparison function to sort triangles based on their score
bool compareTrianglesByScore(const std::pair<Triangle, double>& a, const std::pair<Triangle, double>& b) {
    return a.second < b.second;  // Compare based on scores
}

void printBacktrackCandidateStacks(const std::vector<std::stack<std::vector<Triangle>>>& backtrackCandidateStacks) {
    std::cout << "Printing Backtrack Candidate Stacks:\n";
    for (size_t i = 0; i < backtrackCandidateStacks.size(); ++i) {
        std::cout << "Component " << i << ":\n";

        // Create a temporary stack to hold items while we print and then to restore them back.
        std::stack<std::vector<Triangle>> tempStack;

        std::stack<std::vector<Triangle>> componentStack = backtrackCandidateStacks[i]; // Copy the stack to not affect the original

        int level = 0;
        while (!componentStack.empty()) {
            std::cout << "  Level " << level << ":\n";
            std::vector<Triangle> triangles = componentStack.top();
            componentStack.pop();
            tempStack.push(triangles); // Push into temporary stack

            for (const Triangle& triangle : triangles) {
                std::cout << "    " << triangle.toString() << "\n";
            }

            level++;
        }

        // Restore the items back to componentStack if needed
        while (!tempStack.empty()) {
            componentStack.push(tempStack.top());
            tempStack.pop();
        }
    }
}

AlgorithmState initializeAlgorithmState(
    const std::vector<Component>& tempComponents,
    const std::unordered_set<Triangle>& remainingTriangles,
    const std::unordered_set<Triangle>& usedTriangles) {
    
    AlgorithmState currentState;
    // Existing variable declarations
    currentState.tempComponents = tempComponents;
    currentState.remainingTriangles = remainingTriangles;
    currentState.usedTriangles = usedTriangles;
    currentState.repetitionCounts = std::vector<int>(tempComponents.size(), 1);  // Initialize to 1
    currentState.componentWeights = std::vector<double>(tempComponents.size(), 1.0);  // Initialize to 1.0
    currentState.globalVolumeChangeHistory = std::vector<double>();
    currentState.adjustedCandidateTriangles = std::vector<std::vector<std::pair<double, Triangle>>>(tempComponents.size());
    currentState.backtrackCandidateStacks = std::vector<std::stack<std::vector<Triangle>>>(tempComponents.size());
    currentState.componentsStalled = std::vector<bool>(tempComponents.size(), false);
    currentState.componentIndexStack = std::stack<int>();
    currentState.issueTrianglesStacks = std::vector<std::stack<Triangle>>(tempComponents.size());
    currentState.anyComponentGrown = false;
    currentState.allComponentsStalled = false;
    currentState.triedCombinations = std::unordered_set<std::string>();
    currentState.trianglesThatLedToThisState = std::map<int, std::deque<Triangle>>();
    currentState.candidateTriangles = std::vector<std::vector<std::pair<Triangle, double>>>(tempComponents.size());

    return currentState;

}

enum BacktrackingResult {
    SUCCESS,
    FAILURE,
    EMPTY_STACK,
    BACKTRACKED_PREVIOUSLY
};


BacktrackingResult handleBacktracking(
    AlgorithmState& currentState,
    std::stack<AlgorithmState>& globalBacktrackStack,
    bool& alreadyBacktracked) {

    std::cout << "-----BACKTRACKING NOW-----" << std::endl;
    BacktrackingState backtrackingState;
    // Restore to prevState
    std::string key = generateKeyForCombination(currentState.trianglesThatLedToThisState); // Generate a unique key for this combination
    std::cout << "Key: " << key << std::endl;
    currentState.triedCombinations.insert(key); // Mark this combination as a bad one

    auto tempTriedCombinations = currentState.triedCombinations;
    //originally top, pop, check
    globalBacktrackStack.pop();
    if (globalBacktrackStack.empty()) {
        std::cout << "The global backtrack stack is empty. Cannot proceed." << std::endl;
        return EMPTY_STACK;
    }
    currentState = globalBacktrackStack.top();
    // currentState = globalBacktrackStack.top();
    currentState.triedCombinations = tempTriedCombinations;
    if (!alreadyBacktracked) {
        // printBacktrackCandidateStacks(currentState.backtrackCandidateStacks);
        if (!recursiveBackTracking(globalBacktrackStack, currentState, backtrackingState)) {
            std::cout << "No solution exists!" << std::endl;
            return FAILURE;
        }
        else {
            return SUCCESS;
        }
        alreadyBacktracked = true;
    }
    else {
        return BACKTRACKED_PREVIOUSLY;
    }
}


//Does the seed grow with the temp component? like is the seed attached?
//Why are triangles shrinking so rapidly?
std::pair<bool, std::vector<Component>> growComponentsSynchronously(
    std::vector<Component>& tempComponents,
    std::unordered_set<Triangle>& remainingTriangles,
    TriangleAdjacency& adjacency,
    std::unordered_map<std::string, int>& typeCounts,
    std::unordered_set<Triangle>& usedTriangles,
    bool& dontReseed) {

    std::cout << "Initial number of components: " << tempComponents.size() << std::endl;
    
    AlgorithmState currentState = initializeAlgorithmState(tempComponents, remainingTriangles, usedTriangles);
    std::stack<AlgorithmState> globalBacktrackStack;
    int backtrackCount = 0;
    //place holder for now, it happens when i backtrack once and then backtrack again. I'm going to my very old state instead 
    //of going to my backtrack state which happened most recently after the previous growth didn't work.
    bool alreadyBacktracked = false;


    for (size_t i = 0; i < currentState.tempComponents.size(); ++i) {
        currentState.repetitionCounts[i] = 1;
        currentState.componentWeights[i] = 1.0;
    }

    if (currentState.tempComponents.empty() || currentState.remainingTriangles.empty()) {
        std::cout << "temp components empty or remaining triangles empty" << std::endl;
        exit(1);
        return std::make_pair(currentState.anyComponentGrown, currentState.tempComponents);
    }

    while (true) {
        auto startSection = std::chrono::high_resolution_clock::now();

        std::cout << "----------While loop reset for next iteration-----------" << std::endl;
        std::vector<Triangle> bestTrianglesForAllComponents;  // Initialize here

        // New Code: Reset flags for each iteration
        currentState.anyComponentGrown = false;
        currentState.allComponentsStalled = false;
        
        std::fill(currentState.componentsStalled.begin(), currentState.componentsStalled.end(), false);
        // Retrieve primary candidate triangles
        auto startCandidates = std::chrono::high_resolution_clock::now();
        currentState.candidateTriangles = getCandidateTrianglesForEachComponent(currentState);
        auto stopCandidates = std::chrono::high_resolution_clock::now();
        auto durationCandidates = std::chrono::duration_cast<std::chrono::microseconds>(stopCandidates - startCandidates);
        std::cout << "Candidate triangles generation time: " << durationCandidates.count() << " microseconds" << std::endl;
        printCandidateTrianglesForComponents(currentState.candidateTriangles);
        bool atLeastOneComponentHasCandidates = false;  // Will be set to true if any component has candidate triangles
        // First, check if at least one component has candidate triangles
        for (const auto& candidates : currentState.candidateTriangles) {
            if (!candidates.empty()) {
                atLeastOneComponentHasCandidates = true;
                break;
            }
        }

        // Now, proceed to flag components with no candidate triangles as "stalled", but only if at least one component does have candidate triangles
        if (atLeastOneComponentHasCandidates) {
            bool allComponentsHaveCandidates = true; // Initialize to true
            for (size_t i = 0; i < currentState.candidateTriangles.size(); ++i) {
                if (currentState.candidateTriangles[i].empty()) {
                    std::cout << "Component " << i << " is empty and will be flagged as stalled." << std::endl;
                    currentState.componentIndexStack.push(i);
                    currentState.componentsStalled[i] = true;
                    allComponentsHaveCandidates = false; // Set to false if any component is empty
                } else {
                    std::cout << "Component " << i << " is not empty." << std::endl;
                }
            }

        } else {
            std::cout << "All components have empty candidateTriangles. Exiting..." << std::endl;
            dontReseed = true;
            remainingTriangles = currentState.remainingTriangles;
            usedTriangles = currentState.usedTriangles;
            return std::make_pair(currentState.anyComponentGrown, currentState.tempComponents);
        }


        // Count how many components are stalled
        int numComponentsStalled = std::count(currentState.componentsStalled.begin(), currentState.componentsStalled.end(), true);
        
        //need to look over this
        if (numComponentsStalled > 0) {
            if (globalBacktrackStack.empty()) {
                std::cout << "No alternative paths. Re-seeding..." << std::endl;
                // exit(1);
                // for(auto& tempComponent: currentState.tempComponents) {
                //     std::cout << "Component contains the following triangles:" << std::endl;
                //     for (const auto& triangle : tempComponent.triangles) {
                //         std::cout << "Triangle: " << triangle.toString() << std::endl;
                //     }
                // }
                remainingTriangles = currentState.remainingTriangles;
                usedTriangles = currentState.usedTriangles;
                return std::make_pair(currentState.anyComponentGrown, currentState.tempComponents);
            }
        }

        //If at least one component is stalled, proceed with backtracking
        //ALWAYS POP THEN RETURN. THIS IS NECESSARY BECAUSE WE DONT WANT BAD TRIANGLE GROWTH IN OUR OUTPUT
        if (numComponentsStalled > 0 && !globalBacktrackStack.empty()) {
            BacktrackingResult result = handleBacktracking(currentState, globalBacktrackStack, alreadyBacktracked);
            
            if (result == EMPTY_STACK || result == FAILURE || result == BACKTRACKED_PREVIOUSLY) {
                remainingTriangles = currentState.remainingTriangles;
                usedTriangles = currentState.usedTriangles;
                // Handle empty stack situation
                return std::make_pair(currentState.anyComponentGrown, currentState.tempComponents);
            } else if (result == SUCCESS) {
                // Continue with the next iteration
                continue;
            }
            // Increment the next component to backtrack for the next time.
            // backtrackingState.nextComponentToBacktrack = (backtrackingState.nextComponentToBacktrack + 1) % currentState.tempComponents.size();
        }
        

        std::cout << "Candidate triangles for components size: " << currentState.candidateTriangles.size() << std::endl;

        std::vector<std::vector<std::pair<Triangle, double>>> fullCandidateList = currentState.candidateTriangles;
        auto startSorting = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < currentState.tempComponents.size(); ++i) {
            std::vector<Triangle> allTrianglesForThisComponent;

            // Sort the candidate list for this component based on scores
            std::sort(fullCandidateList[i].begin(), fullCandidateList[i].end(), compareTrianglesByScore);

            // Extract just the triangles from the sorted candidate list
            for (const auto& pair : fullCandidateList[i]) {
                allTrianglesForThisComponent.push_back(pair.first);
            }

            // Push the sorted list onto the stack for this component
            if (!allTrianglesForThisComponent.empty()) {
                currentState.backtrackCandidateStacks[i].push(allTrianglesForThisComponent);
                std::cout << "Pushed " << allTrianglesForThisComponent.size() << " triangles onto backtrack stack for component " << i << std::endl;
            }
        }
        auto stopSorting = std::chrono::high_resolution_clock::now();
        auto durationSorting = std::chrono::duration_cast<std::chrono::microseconds>(stopSorting - startSorting);
        std::cout << "Triangle sorting time: " << durationSorting.count() << " microseconds" << std::endl;

        std::vector<double> minConvexHullChanges(currentState.tempComponents.size(), std::numeric_limits<double>::infinity());
        for (auto& vec : currentState.adjustedCandidateTriangles) {
            vec.clear();
        }
        
        std::vector<int> hypotheticalNewZi(currentState.tempComponents.size());
        std::vector<double> hypotheticalNewWi(currentState.tempComponents.size());
        for (size_t j = 0; j < currentState.tempComponents.size(); ++j) {
            int oldZi = currentState.repetitionCounts[j];
            hypotheticalNewZi[j] = oldZi + 1;
            double oldWi = currentState.componentWeights[j];
            hypotheticalNewWi[j] = calculateWeight(hypotheticalNewZi[j]);
        }

        std::vector<std::vector<std::pair<double, Triangle>>> bestTrianglesForComponents(currentState.tempComponents.size());

        for(size_t i = 0; i < currentState.tempComponents.size(); ++i) {
            int repetitionCount = currentState.repetitionCounts[i]; // Z_i in Eqn. (4)
            double weightImpact = currentState.componentWeights[i]; // w_i in Eqn. (4)
            
            std::cout << "Component " << i << " repetition count : " << repetitionCount << std::endl;
            std::cout << "Component " << i << " weight impact: " << weightImpact << std::endl;

            for(const auto& candidatePair : currentState.candidateTriangles[i]) {
                const Triangle& candidate = candidatePair.first;
                
                // Check for contiguity before even calculating convex hull changes
                if (!isValidAfterAddingTriangle(currentState.tempComponents[i], candidate)) {
                    continue;  // Skip non-contiguous triangles
                }
                
                double change = calculateConvexHullChange(currentState.tempComponents[i], candidate);
                change *= hypotheticalNewWi[i];  // Use the updated weight
                
                // Add to bestTriangles with the change for sorting later
                bestTrianglesForComponents[i].emplace_back(change, candidate);
            }

            // Sort the triangles within each component based on the minimizing 'change'
            std::sort(bestTrianglesForComponents[i].begin(), bestTrianglesForComponents[i].end(), 
                    [](const std::pair<double, Triangle>& a, const std::pair<double, Triangle>& b) {
                        return a.first < b.first;
                    });
        }


        std::cout << "ABOUT TO SYNCHRONIZE GROWTH" << std::endl;
        
        // auto bestTrianglesForComponents = synchronizeGrowth(currentState, minConvexHullChanges, true);
        // printBestTrianglesForComponents(bestTrianglesForComponents);
        bool shouldBreakOuterLoop = false;
        for (size_t i = 0; i < currentState.tempComponents.size(); ++i) {
            // Validate that all best triangles can be added and are not already used
            bool allUpdatesAreValid = true;
            bool allBestTrianglesAreValid = true;
            for (size_t j = 0; j < currentState.tempComponents.size(); ++j) {
                // Check if bestTrianglesForComponents[j] is empty
                Triangle bestTriangleForValidation = bestTrianglesForComponents[j][0].second;
                bestTrianglesForAllComponents.push_back(bestTriangleForValidation);  // Populate here
                
                if (!isValidAfterAddingTriangle(currentState.tempComponents[j], bestTriangleForValidation) || 
                    !isValidWeight(hypotheticalNewWi[j], currentState.componentWeights[j])) {
                    allUpdatesAreValid = false;
                    break;
                }

                if (currentState.usedTriangles.find(bestTriangleForValidation) != currentState.usedTriangles.end()) {
                    std::cout << "Triangle " << bestTriangleForValidation.toString() << " is already used." << std::endl;
                    allBestTrianglesAreValid = false;
                    break;
                }
            }
            if (allUpdatesAreValid && allBestTrianglesAreValid) {
                for (size_t j = 0; j < currentState.tempComponents.size(); ++j) {
                    Triangle bestTriangle = bestTrianglesForComponents[j][0].second;
                    std::cout << "Adding triangle " << bestTriangle.toString() << " to component " << j << std::endl;
                    currentState.tempComponents[j].addTriangle(bestTriangle);
                    currentState.repetitionCounts[j] = hypotheticalNewZi[j];
                    currentState.componentWeights[j] = hypotheticalNewWi[j];
                    double convexHullChange = minConvexHullChanges[j];
                    currentState.anyComponentGrown = true;
                    currentState.usedTriangles.insert(bestTriangle);
                    currentState.remainingTriangles.erase(bestTriangle);
                    currentState.trianglesThatLedToThisState[j].push_back(bestTriangle);
                    currentState.globalVolumeChangeHistory.emplace_back(convexHullChange);
                    bestTrianglesForComponents[j].erase(bestTrianglesForComponents[j].begin());
                    shouldBreakOuterLoop = true;  // Set the flag to true
                }
                globalBacktrackStack.push(currentState);
            }
            std::cout << "Temp component " << i << " triangles after trying to grow: " << currentState.tempComponents[i].triangles.size() << std::endl;
            if (shouldBreakOuterLoop) {
                break;  // This will break out of the outer loop
            }
            
        }
        auto stopSection = std::chrono::high_resolution_clock::now();
        auto durationSection = std::chrono::duration_cast<std::chrono::microseconds>(stopSection - startSection);
        std::cout << "Total time for this section: " << durationSection.count() << " microseconds" << std::endl;
    }
}


bool validateFinalState(const std::vector<Component>& components, const Cluster& cluster) {
    std::unordered_set<Triangle> allTrianglesInComponents;

    for (const auto& component : components) {
        // Check Property i: Ensure triangles in each component are contiguous
        if (!areTrianglesContiguous(component.triangles)) {
            std::cout << "Warning: Component is not contiguous!" << std::endl;
            return false;
        }

        // Add triangles to the set for later checks
        allTrianglesInComponents.insert(component.triangles.begin(), component.triangles.end());
    }

    // std::cout << "Total triangles in all components: " << allTrianglesInComponents.size() << std::endl;
    // std::cout << "Total triangles in cluster: " << cluster.triangles.size() << std::endl;

    // Check Property ii: Ensure all triangles are accounted for
    if (allTrianglesInComponents.size() != cluster.triangles.size()) {
        std::cout << "Warning: Some triangles are missing!" << std::endl;

        // Optionally, print the specific triangles that are missing
        for (const auto& triangle : cluster.triangles) {
            if (allTrianglesInComponents.find(triangle) == allTrianglesInComponents.end()) {
                std::cout << "Missing Triangle: " << triangle.toString() << std::endl;
            }
        }
        return false;
    }

    return true;
}

bool isComponentInList(const Component& component, const std::vector<Component>& componentList) {
    for (const auto& existingComponent : componentList) {
        if (component == existingComponent) {  // Replace '==' with your equality check
            return true;
        }
    }
    return false;
}

std::vector<Component> randomizedGrowthOptimization(const Cluster& cluster,
                                                    TriangleAdjacency& adjacency,
                                                    std::vector<double> weights,
                                                    double tau_S,
                                                    int& iteration) {
    std::cout << "Entering randomizedGrowthOptimization function" << std::endl;

    std::unordered_map<ComponentType, double> componentWeights;  // To store the weights of each component type

    ComponentType deepestType;  // Declare it here or outside the loop based on your need
    int deepestTypeDepth = -1;  // Initialize to -1 as a sentinel value
    std::vector<Component> finalComponents;
    std::unordered_set<Triangle> remainingTriangles(cluster.triangles.begin(), cluster.triangles.end());
    // std::unordered_map<std::string, int> typePriority;  // For re-seeding priority based on types
    std::unordered_map<std::string, int> typeCounts;
    std::unordered_set<Triangle> usedTriangles;  // <-- Declare here to keep track of used triangles
    static bool dontReseed = false;

    std::cout << "-----iteration: " << iteration << " starting----- " << std::endl;
    std::cout << "Number of triangles in cluster: " << cluster.triangles.size() << std::endl;
    bool needsReseeding = false;
    std::vector<Component> tempComponents;
    std::vector<Component> deepestComponents;
    std::vector<Triangle> seeds;
    // Regular seeding and initialization
    seeds = findSeeds(remainingTriangles, seeds, weights, tau_S);
    if (seeds.empty()) {
        std::cout << "No seeds found. Exiting." << std::endl;
        exit(1);
    }
    tempComponents = initializeTempComponents(seeds, remainingTriangles);
    if(remainingTriangles.empty()) {
        for (const auto& tempComponent : tempComponents) {
                finalComponents.push_back(tempComponent);
        }
    }
    while (!remainingTriangles.empty()) {
        std::cout << "Triangle count before growing in remainingTriangles: " << remainingTriangles.size() << std::endl;
        std::cout << "Triangle count in seeds: " << seeds.size() << std::endl;
        std::pair<bool, std::vector<Component>> result = growComponentsSynchronously(
            tempComponents, remainingTriangles, adjacency,
            typeCounts, usedTriangles, dontReseed
        );
        bool canGrow = result.first;
        deepestComponents = result.second;
        // If we couldn't grow, use the deepest components found so far
        // and then re-seed with remaining types.
        finalComponents.insert(finalComponents.end(), deepestComponents.begin(), deepestComponents.end());
        // Reset tempComponents for the next iteration
        std::cout << "temp components size: " << tempComponents.size() << std::endl;
        tempComponents.clear();
        usedTriangles.clear();
        seeds.clear();
        needsReseeding = true;
        // logic to re-seed with remaining types (e.g., C's if AB was the deepest)
        seeds = findSeeds(remainingTriangles, seeds, weights, tau_S);

        // std::cout << "seed count after reseeding: " << seeds.size() << std::endl;
        std::vector<Component> newTempComponents = initializeTempComponents(seeds, remainingTriangles);
        // std::cout << "Triangle count after reseeding in remainingTriangles: " << remainingTriangles.size() << std::endl;
        tempComponents.insert(tempComponents.end(), newTempComponents.begin(), newTempComponents.end());
        // std::cout << "RESEEDED PROPERLY" << std::endl;
        if (remainingTriangles.empty()) {
            for (const auto& tempComponent : tempComponents) {
                // if (!isComponentInList(tempComponent, finalComponents)) {
                    finalComponents.push_back(tempComponent);
                // }
            }
            // std::cout << "Added unique remaining tempComponents to finalComponents as remainingTriangles is empty." << std::endl;
            break;  // Exit the while loop
        }
        if(tempComponents.empty()) {
            std::cout << "tempComponents is empty. Exiting." << std::endl;
            break;
        }
    }

    // Before returning finalComponents, perform the merging operation
    std::cout << "Before merging: " << finalComponents.size() << std::endl;
    for (const auto& component : finalComponents) {
        std::cout << "Component has " << component.triangles.size() << " triangles." << std::endl;
        for(const auto& triangle : component.triangles) {
            std::cout << "Triangle " << triangle.toString() << " is in component." << std::endl;
        }
    }
    mergeOverlappingComponents(finalComponents, adjacency);
    std::cout << "Triangle count at the start: " << cluster.triangles.size() << std::endl;
    std::cout << "After merging: " << finalComponents.size() << std::endl;

    // Final validation
    if (!validateFinalState(finalComponents, cluster)) {
        std::cout << "Validation failed." << std::endl;
        exit(1);  // Exits the program with a status of 1
    }

    seeds.clear();

    return finalComponents;
}

void optimizeAndCategorizeComponents(std::vector<Cluster>& allClusters,
                                     TriangleAdjacency& adjacency,
                                     std::vector<double> weights,
                                     double tau_S,
                                     std::vector<Component>& allFinalComponents) {  // Passed by reference
    std::unordered_set<ComponentType> C;  // Global set of unique component types
    std::unordered_map<ComponentType, std::vector<Component>> Z;  // Global instances of each component type
    std::unordered_map<ComponentType, double> componentWeights;  // Global weights of each component type
    int iteration = 0;
    int componentID = 0; // ID counter for components

    // This loop represents multiple iterations of randomized growth optimization
    for (Cluster& cluster : allClusters) {
        std::vector<Component> finalComponentsForCluster = randomizedGrowthOptimization(
            cluster, adjacency, weights, tau_S, iteration
        );

        allFinalComponents.insert(
            allFinalComponents.end(),
            finalComponentsForCluster.begin(),
            finalComponentsForCluster.end()
        );
    }
    std::cout << "ALL DONE WITH GROWTH FUNCTION" << std::endl;

        // Now update C and Z based on finalComponentsForCluster
        // for (const Component& component : finalComponentsForCluster) {
        //     ComponentType type = identifyComponentType(component);
        //     if (C.find(type) != C.end()) {
        //         Z[type].push_back(component);
        //         componentWeights[type] = 1.0 / std::pow(Z[type].size(), 2);
        //     } else {
        //         C.insert(type);
        //         Z[type] = {component};
        //         componentWeights[type] = 1.0;
        //     }
        // }
    // }
}


// Function to assign contrasting hues to clusters
std::vector<float> assignContrastingHues(int numClusters) {
    std::vector<float> hues;
    float currentHue = 0.0f;
    float hueStep = 360.0f / numClusters;

    // Assign every other hue first to maximize difference between adjacent clusters
    for (int i = 0; i < numClusters; i += 2) {
        hues.push_back(currentHue);
        currentHue += hueStep * 2;
    }

    // Now fill in the gaps
    currentHue = hueStep;
    for (int i = 1; i < numClusters; i += 2) {
        hues.push_back(currentHue);
        currentHue += hueStep * 2;
    }

    // Shuffle to randomize
    std::random_shuffle(hues.begin(), hues.end());

    return hues;
}


Vector3D HSVtoRGB(float H, float S, float V) {
    if (S == 0) {
        return Vector3D(V, V, V);  // It's a shade of gray
    }

    H /= 60.0f;  // Split the hue into sextant and fractional parts
    int i = static_cast<int>(H);
    float F = H - i;

    float P = V * (1 - S);
    float Q = V * (1 - S * F);
    float T = V * (1 - S * (1 - F));

    Vector3D color;

    switch (i) {
        case 0:
            color = Vector3D(V, T, P);
            break;
        case 1:
            color = Vector3D(Q, V, P);
            break;
        case 2:
            color = Vector3D(P, V, T);
            break;
        case 3:
            color = Vector3D(P, Q, V);
            break;
        case 4:
            color = Vector3D(T, P, V);
            break;
        default:
            color = Vector3D(V, P, Q);
            break;
    }

    return color;
}

// Function to assign contrasting hues, saturations, and values to clusters
void assignContrastingColors(std::vector<Cluster>& clusters) {
    int numGroups = std::sqrt(clusters.size());
    float hueStep = 360.0f / numGroups;

    std::vector<float> saturations = {1.0f, 0.8f, 0.6f};
    std::vector<float> values = {1.0f, 0.8f};

    int groupIndex = 0;
    int saturationIndex = 0;
    int valueIndex = 0;

    for (int i = 0; i < clusters.size(); ++i) {
        if (i % numGroups == 0) {
            // Change the group and update saturation and value
            groupIndex++;
            saturationIndex = (saturationIndex + 1) % saturations.size();
            valueIndex = (valueIndex + 1) % values.size();
        }

        float currentHue = static_cast<float>(groupIndex) * hueStep;
        float currentSaturation = saturations[saturationIndex];
        float currentValue = values[valueIndex];

        Vector3D assignedColor = HSVtoRGB(currentHue, currentSaturation, currentValue);

        for (auto& triangle : clusters[i].triangles) {
            triangle.setColor(assignedColor);
        }
    }
}

void assignContrastingColorsToComponents(std::vector<Component>& components) {
    int numGroups = std::sqrt(components.size());
    float hueStep = 360.0f / numGroups;

    std::vector<float> saturations = {1.0f, 0.8f, 0.6f};
    std::vector<float> values = {1.0f, 0.8f};

    int groupIndex = 0;
    int saturationIndex = 0;
    int valueIndex = 0;

    for (int i = 0; i < components.size(); ++i) {
        if (i % numGroups == 0) {
            // Change the group and update saturation and value
            groupIndex++;
            saturationIndex = (saturationIndex + 1) % saturations.size();
            valueIndex = (valueIndex + 1) % values.size();
        }

        float currentHue = static_cast<float>(groupIndex) * hueStep;
        float currentSaturation = saturations[saturationIndex];
        float currentValue = values[valueIndex];

        Vector3D assignedColor = HSVtoRGB(currentHue, currentSaturation, currentValue);

        for (auto& triangle : components[i].triangles) {
            triangle.setColor(assignedColor);
        }
    }
}



void validateTotalTriangles(const std::vector<Cluster>& clusters, size_t initialCount) {
    size_t currentCount = 0;
    for (const auto& cluster : clusters) {
        currentCount += cluster.triangles.size();
    }
    if (currentCount != initialCount) {
        std::cout << "Triangle count mismatch. Initial: " << initialCount << ", Current: " << currentCount << std::endl;
    } else {
        std::cout << "Triangle count validated. Total: " << currentCount << std::endl;
    }
}

void checkForDuplicateTriangles(const Cluster& cluster) {
    std::unordered_map<Triangle, int> triangleCount;
    for (const auto& triangle : cluster.triangles) {
        triangleCount[triangle]++;
    }

    for (const auto& [triangle, count] : triangleCount) {
        if (count > 1) {
            std::cout << "Duplicate triangle detected: " << triangle.toString() << " appears " << count << " times." << std::endl;
        }
    }
}
void checkForDuplicates(const std::vector<Cluster>& clusters) {
    // Check for duplicates within each cluster
    for (const auto& cluster : clusters) {
        std::unordered_set<Triangle> uniqueTriangles(cluster.triangles.begin(), cluster.triangles.end());
        if (uniqueTriangles.size() != cluster.triangles.size()) {
            std::cerr << "Duplicate triangles found in a cluster" << std::endl;
        }
    }

    // Check for duplicates across clusters
    std::unordered_set<Triangle> allTriangles;
    for (const auto& cluster : clusters) {
        for (const auto& triangle : cluster.triangles) {
            if (allTriangles.count(triangle) > 0) {
                std::cerr << "Duplicate triangle detected across clusters: " << triangle.toString() << std::endl;
            }
            allTriangles.insert(triangle);
        }
    }
}


void initializeOriginalClusterOfTriangle(std::unordered_map<Triangle, size_t>& originalClusterOfTriangle, const std::vector<Cluster>& clusters) {
    originalClusterOfTriangle.clear();  // Clear the existing map
    for (size_t i = 0; i < clusters.size(); ++i) {
        for (const auto& triangle : clusters[i].triangles) {
            originalClusterOfTriangle[triangle] = i;
        }
    }
}

void refineClusters(std::vector<Cluster>& clusters, double tau_N, double spatialHashCellSize) {
    int iteration = 0;
    SpatialHash spatialHash(spatialHashCellSize);
    auto start = std::chrono::high_resolution_clock::now();
    spatialHash.precomputeTriangleHashes(clusters);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Time taken for precomputeTriangleHashes: " << elapsed.count() << " seconds" << std::endl;
    // 1. Initialize the TriangleAdjacency object
    TriangleAdjacency triangleAdjacency;
    for (const auto& cluster : clusters) {
        for (const auto& triangle : cluster.triangles) {
            triangleAdjacency.addTriangle(triangle);
        }
    }
    std::unordered_map<Triangle, std::vector<Triangle>> precomputedNeighbors;
    auto potentialNeighborsTimeStart = std::chrono::high_resolution_clock::now();
    for (const auto& cluster : clusters) {
        for (const auto& triangle : cluster.triangles) {
            precomputedNeighbors[triangle] = spatialHash.getPotentialNeighbors(triangle, triangleAdjacency);
        }
    }
    auto potentialNeighborsTimeEnd = std::chrono::high_resolution_clock::now();

    auto potentialNeighborsDuration = std::chrono::duration_cast<std::chrono::milliseconds>(potentialNeighborsTimeEnd - potentialNeighborsTimeStart);
    std::cout << "Time for finding potential neighbors: " << potentialNeighborsDuration.count() << " ms\n";

    while (true) {

        std::cout << "Starting iteration " << iteration << " with " << clusters.size() << " clusters." << std::endl;

        int mergesInThisIteration = 0;
        int splitsInThisIteration = 0;
        size_t pairsProcessed = 0;  // Initialize counter for this iteration
        std::unordered_set<Triangle> alreadyAddedToMaps;  // New line


        // Declare a variable to accumulate the total time for the specific loop
        auto totalStartTime = std::chrono::high_resolution_clock::now();
        std::chrono::milliseconds totalLoopTime(0);  // Initialize to zero


        std::map<std::pair<size_t, size_t>, std::vector<Triangle>> mergeMap;
        std::map<std::pair<size_t, size_t>, std::vector<Triangle>> splitMap;
        std::unordered_map<int, std::unordered_set<int>> neighboringClustersCache;

        for (size_t i = 0; i < clusters.size(); ++i) {
            if (clusters[i].triangles.empty()) continue;

            if (neighboringClustersCache.find(i) == neighboringClustersCache.end()) {
                neighboringClustersCache[i] = spatialHash.getNeighboringClustersForCluster(clusters[i], triangleAdjacency);
            }

            auto& neighboringClusterIndices = neighboringClustersCache[i];

            std::vector<int> neighboringClusterIndicesVec(neighboringClusterIndices.begin(), neighboringClusterIndices.end());
            for (size_t idx = 0; idx < neighboringClusterIndicesVec.size(); ++idx) {
                int j = neighboringClusterIndicesVec[idx];
                if (i == j || clusters[j].triangles.empty()) continue;

                std::unordered_set<Triangle> localMergeList;
                std::unordered_set<Triangle> tbTriangles(clusters[j].triangles.begin(), clusters[j].triangles.end());

                double similarity = ATaTb(clusters[i], clusters[j], precomputedNeighbors, triangleAdjacency, tbTriangles, tau_N);
                std::unordered_set<Triangle> toBeMerged;
                std::unordered_set<Triangle> toBeSplit;

                if (similarity >= tau_N) {
                    for (const auto& triangleA : clusters[i].triangles) {
                        bool foundAdjacent = false;
                        auto it = precomputedNeighbors.find(triangleA);
                        if (it != precomputedNeighbors.end()) {
                            auto& cachedPotentialNeighbors = it->second;
                            for (const auto& triangleB : cachedPotentialNeighbors) {
                                if (triangleA.isAdjacent(triangleB, triangleAdjacency)) {
                                    // Both triangles meet the merge condition, so add them to the list
                                    toBeMerged.insert(triangleA);
                                    toBeMerged.insert(triangleB);  // Insert triangleB as well
                                    foundAdjacent = true;
                                    // break;
                                }
                            }
                            if (!foundAdjacent) {
                                toBeSplit.insert(triangleA);
                            }
                        }
                    }
                    for (const auto& triangleB : clusters[j].triangles) {
                            if (toBeMerged.count(triangleB) == 0) {
                                toBeSplit.insert(triangleB);
                            }
                    }
                    if (!toBeMerged.empty()) {
                        for (const auto& triangle : toBeMerged) {
                            if (alreadyAddedToMaps.count(triangle) == 0) {
                                mergeMap[{i, j}].emplace_back(triangle);
                                alreadyAddedToMaps.insert(triangle);
                            }
                        }
                    }
                    
                    if (!toBeSplit.empty()) {
                        for (const auto& triangle : toBeSplit) {
                            if (alreadyAddedToMaps.count(triangle) == 0) {
                                splitMap[{i, j}].emplace_back(triangle);
                                alreadyAddedToMaps.insert(triangle);
                            }
                        }
                    }
                }
            }
        }
        auto totalEndTime = std::chrono::high_resolution_clock::now();  // End time measurement
        auto totalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(totalEndTime - totalStartTime);

        std::cout << "Total time for entire program: " << totalDuration.count() << " ms\n";

        size_t initialCount = 0;
        for (const auto& cluster : clusters) {
            initialCount += cluster.triangles.size();
        }
        std::cout << "Initial Triangle Count: " << initialCount << std::endl;

        // Initialize a variable to keep track of the total number of triangles after all merges
        std::unordered_set<Triangle> mergedTriangles;  // To keep track of merged triangles
        std::unordered_set<size_t> mergedClusters;
        std::unordered_map<Triangle, size_t> originalClusterOfTriangle;  // Triangle -> Original cluster index

        initializeOriginalClusterOfTriangle(originalClusterOfTriangle, clusters);

        // Create a temporary copy of the clusters
        std::vector<Cluster> tempClusters = clusters;
        clusters.clear();

        std::cout << "Size of mergeMap before first merge: " << mergeMap.size() << std::endl;
        std::cout << "Size of splitMap before first merge: " << splitMap.size() << std::endl;


        for (auto& entry : mergeMap) {
            auto& keyPair = entry.first;
            auto& trianglesToMerge = entry.second;
            size_t sourceClusterIndex = keyPair.first;
            size_t targetClusterIndex = keyPair.second;

            auto& sourceTriangles = tempClusters[sourceClusterIndex].triangles;
            auto& targetTriangles = tempClusters[targetClusterIndex].triangles;
            std::unordered_set<Triangle> targetClusterTriangles(targetTriangles.begin(), targetTriangles.end());

            // Merge logic
            for (const auto& triangle : trianglesToMerge) {
                // Check if this triangle has already been merged
                if (mergedTriangles.find(triangle) != mergedTriangles.end()) {
                    continue;  // Skip this triangle
                }

                // Check if the triangle exists in the source cluster before removing it
                if (std::find(sourceTriangles.begin(), sourceTriangles.end(), triangle) != sourceTriangles.end()) {
                    // Remove the triangle from the source cluster
                    auto it = std::remove(sourceTriangles.begin(), sourceTriangles.end(), triangle);
                    sourceTriangles.erase(it, sourceTriangles.end());

                    // Add the triangle to the target cluster if it doesn't already exist there
                    if (targetClusterTriangles.find(triangle) == targetClusterTriangles.end()) {
                        targetTriangles.emplace_back(triangle);
                        targetClusterTriangles.emplace(triangle);
                    }

                    // Mark this triangle as merged
                    mergedTriangles.insert(triangle);
                }
            }
            mergedClusters.insert(targetClusterIndex);
        }

        clusters = tempClusters;  // Update clusters with merged state

        // After all merge operations are complete and clusters have been updated
        for (size_t i = 0; i < clusters.size(); ++i) {
            for (const auto& triangle : clusters[i].triangles) {
                originalClusterOfTriangle[triangle] = i;
            }
        }

        std::unordered_set<Triangle> movedTriangles;  // To keep track of triangles that have already been moved

        std::unordered_set<Triangle> splitTriangles;  // To keep track of triangles that have been split
        // std::chrono::milliseconds totalLoopSplitTime(0);
        // auto loopSplitStart = std::chrono::high_resolution_clock::now();  // Start time measurement
        std :: cout << "cluster size before split: " << clusters.size() << std::endl;

        size_t initialTriangleCount = 0;
        for (const auto& cluster : clusters) {
            initialTriangleCount += cluster.triangles.size();
        }

        validateTotalTriangles(clusters, initialTriangleCount);

        // Remove empty clusters after merges and splits
        auto mergeClear = std::remove_if(clusters.begin(), clusters.end(),
            [](const Cluster& cluster) { return cluster.triangles.empty(); });
            
        std::cout << "Number of empty clusters: " << std::distance(mergeClear, clusters.end()) << std::endl;

        clusters.erase(mergeClear, clusters.end());

        checkForDuplicates(clusters);
        size_t trianglesAdded = 0;
        size_t trianglesRemoved = 0;

       // Your existing random generator and other variables remain unchanged
        // Initialize random number generator and other containers
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distr(0, clusters.size() - 1);

        std::unordered_set<size_t> emptiedClusters; // Keep track of emptied clusters

        std::cout << "Initial number of clusters: " << clusters.size() << std::endl;

        std::unordered_set<Triangle> processedTriangles;

        for (auto& entry : splitMap) {
            auto& trianglesToSplit = entry.second;
            size_t originalClusterIndex = originalClusterOfTriangle[trianglesToSplit[0]];

            if (emptiedClusters.count(originalClusterIndex) > 0) {
                continue;
            }

            // Remove triangles from their original cluster
            auto& originalClusterTriangles = clusters[originalClusterIndex].triangles;
            for (const auto& triangle : trianglesToSplit) {
                if (processedTriangles.count(triangle) == 0) {
                    auto it = std::remove(originalClusterTriangles.begin(), originalClusterTriangles.end(), triangle);
                    if (it != originalClusterTriangles.end()) {
                        originalClusterTriangles.erase(it, originalClusterTriangles.end());
                        processedTriangles.insert(triangle);
                    }
                }
            }

            // Check if the original cluster is now empty
            if (originalClusterTriangles.empty()) {
                emptiedClusters.insert(originalClusterIndex);
            }

            // Create a new cluster and add it to clusters
            Cluster newCluster;
            for (const auto& triangle : trianglesToSplit) {
                if (processedTriangles.count(triangle) > 0) {
                    newCluster.triangles.emplace_back(triangle);
                }
            }

            // Only add non-empty new clusters
            if (!newCluster.triangles.empty()) {
                clusters.push_back(newCluster);

                // Update the originalClusterOfTriangle map
                for (const auto& triangle : newCluster.triangles) {
                    originalClusterOfTriangle.emplace(triangle, clusters.size() - 1);
                }
            }
        }



        auto itSplit = std::remove_if(clusters.begin(), clusters.end(),
            [](const Cluster& cluster) { return cluster.triangles.empty(); });

        std::cout << "Number of empty clusters: " << std::distance(itSplit, clusters.end()) << std::endl;

        clusters.erase(itSplit, clusters.end());

        std::cout << "Final number of clusters: " << clusters.size() << std::endl;


        // std::random_device rd;
        // std::mt19937 gen(rd());
        // std::uniform_int_distribution<> distr(0, clusters.size() - 1);
        // std::unordered_set<size_t> emptiedClusters;  // Keep track of emptied clusters
        // std::unordered_map<Triangle, size_t> triangleMovedToCluster;

        // for (auto& entry : splitMap) {
        //     auto& trianglesToSplit = entry.second;
        //     size_t originalClusterIndex = originalClusterOfTriangle[trianglesToSplit[0]];

        //     if (emptiedClusters.count(originalClusterIndex) > 0) {
        //         continue;
        //     }

        //     auto& originalClusterTriangles = clusters[originalClusterIndex].triangles;

        //     size_t notMovedCount = 0;
        //     std::vector<Triangle> successfullyMoved;

        //     for (const auto& triangle : originalClusterTriangles) {
        //         if (movedTriangles.count(triangle) == 0) {
        //             successfullyMoved.push_back(triangle);
        //             movedTriangles.insert(triangle);  // Mark as moved, but don't actually move yet
        //         } else {
        //             notMovedCount++;
        //         }
        //     }

        //     emptiedClusters.insert(originalClusterIndex);  // Mark this cluster as to be emptied, but don't actually clear it yet
        // }


        // // Remove empty clusters after merges and splits
        // auto itSplit = std::remove_if(clusters.begin(), clusters.end(),
        //     [](const Cluster& cluster) { return cluster.triangles.empty(); });
            
        // std::cout << "Number of empty clusters: " << std::distance(itSplit, clusters.end()) << std::endl;

        // clusters.erase(itSplit, clusters.end());


        // Your function to check for duplicates (Assuming you have this implemented)
        // checkForDuplicates(clusters);



        size_t finalTriangleCount = 0;
        for (const auto& cluster : clusters) {
            finalTriangleCount += cluster.triangles.size();
        }

        // Verify
        if (initialTriangleCount != finalTriangleCount || trianglesAdded != trianglesRemoved) {
            std::cerr << "ERROR: Triangle count mismatch. Initial: " << initialTriangleCount << ", Final: " << finalTriangleCount << std::endl;
        }

        validateTotalTriangles(clusters, initialTriangleCount);

        // std::cout << "Total time taken: " << totalLoopSplitTime.count() << " microseconds" << std::endl;

        std::cout << "After " << iteration << " iterations, merged " << mergeMap.size() << " triangles and split " << splitMap.size() << " triangles." << std::endl;

        // After you calculate the merged and split triangles for this iteration
        int mergedAndSplitThisIteration = mergeMap.size() + splitMap.size();
        std::cout << "Merged and split this iteration: " << mergedAndSplitThisIteration << std::endl;

        if (mergeMap.empty() && splitMap.empty()) {
                std::cout << "Convergence achieved. Stopping..." << std::endl;
                break;
        }
        else {
            std::cout << "convergence not achieved based on change rate. Continuing..." << std::endl;
        }

        mergeMap.clear();
        splitMap.clear();
        spatialHash.clear();
        spatialHash.precomputeTriangleHashes(clusters);

        // find out how many clusters there are now
        std::cout << "After erase and spatial clearing: " << iteration << " iterations, there are " << clusters.size() << " clusters." << std::endl;


        // After the merge and split operations
        std::unordered_map<Triangle, int> triangleCount;
        for (const auto& cluster : clusters) {
            for (const auto& triangle : cluster.triangles) {
                triangleCount[triangle]++;
            }
        }

        // Consistency check
        size_t currentTriangleCount = 0;
        for (const auto& cluster : clusters) {
            currentTriangleCount += cluster.triangles.size();
        }

        if (currentTriangleCount != initialTriangleCount) {
            std::cout << "ERROR: Triangle count mismatch. Initial: " << initialTriangleCount << ", Current: " << currentTriangleCount << std::endl;
            break;  // Exiting the loop if there's a mismatch
        }

        iteration++;
    }

}


std::vector<Cluster> createSearchSpaces(const std::vector<Triangle>& triangles, std::vector<double> weightClusters, std::vector<double> weightComponents, 
                                        double tau_S, double spatialHashCellSize) {    // 1. Initial clustering based on shape similarity
    std::cout << "1. Initial clustering based on shape similarity" << std::endl;
    std::vector<double> within_cluster_stitj;
    std::vector<double> between_cluster_stitj;

    // Your clustering function call
    // std::vector<Cluster> clusters = initialClusteringByShapeSimilarity(triangles, weights, tau_S, maxArea, maxEdgeMax, maxEdgeMin, within_cluster_stitj, between_cluster_stitj);

    std::vector<Cluster> clusters = initialClusteringByShapeSimilarity(triangles, weightClusters, tau_S);

    std::cout << "Number of initial clusters: " << clusters.size() << std::endl;
    std::cout << "Initial clustering done. Number of clusters: " << clusters.size() << std::endl;
    
    // 3. Iterative merging and splitting based on adjacency similarity
    std::cout << "3. Iterative merging and splitting based on adjacency similarity" << std::endl;
    double tau_N = 0.5;  // Threshold for adjacency similarity
    refineClusters(clusters, tau_N, spatialHashCellSize);

    std::cout << "Number of clusters after refining: " << clusters.size() << std::endl;

    std::vector<Component> allFinalComponents;

    // // 4. Randomized Growth Optimization
    std::cout << "4. Randomized Growth Optimization" << std::endl;
    TriangleAdjacency adjacencyInstance; // Assume you've created and possibly initialized this instance
    for (const Cluster& cluster : clusters) {
        for (const Triangle& triangle : cluster.triangles) {
            adjacencyInstance.addTriangle(triangle);
        }
    }

    optimizeAndCategorizeComponents(clusters, adjacencyInstance, weightComponents, tau_S, allFinalComponents);
    
    std::cout << "Total components after optimize and categorize components: " << allFinalComponents.size() << std::endl;


    std::cout << "Before assignContrastingColors" << std::endl;
    assignContrastingColors(clusters);
    std::cout << "After assignContrastingColors" << std::endl;

    std::cout << "Before assignContrastingColorsToComponents" << std::endl;
    assignContrastingColorsToComponents(allFinalComponents);
    std::cout << "After assignContrastingColorsToComponents" << std::endl;

    std::cout << "Before saveToOBJ" << std::endl;
    saveToOBJ(clusters, "clusters");
    std::cout << "After saveToOBJ" << std::endl;

    std::cout << "Before saveComponentsToOBJ" << std::endl;
    saveComponentsToOBJ(allFinalComponents, "components");
    std::cout << "After saveComponentsToOBJ" << std::endl;

    std::cout << "Saved refined clusters to OBJ for visualization." << std::endl;

}



// void testSpatialHash(const std::vector<Triangle>& triangles) {

//     SpatialHash spatialHash(0.1);  // Assuming 1 is the cell size you want

//     // Ensure there are at least 3 triangles for the test
//     if (triangles.size() < 3) {
//         std::cout << "Not enough triangles for the test." << std::endl;
//         return;
//     }

//     // Use the first 3 triangles for the test
//     std::cout << "Before triangle 1 assignment" << std::endl;
//     Triangle t1 = triangles[0];
//     std::cout << "After triangle 1 assignment" << std::endl;
//     Triangle t2 = triangles[1];
//     Triangle t3 = triangles[2];

//     spatialHash.insert(t1, 0);
//     spatialHash.insert(t2, 1);
//     spatialHash.insert(t3, 2);

//     auto neighborsT1 = spatialHash.getPotentialNeighbors(t1);
//     auto neighborsT2 = spatialHash.getPotentialNeighbors(t2);
//     auto neighborsT3 = spatialHash.getPotentialNeighbors(t3);

//     std::cout << "Neighbors for T1: " << neighborsT1.size() << std::endl;
//     std::cout << "Neighbors for T2: " << neighborsT2.size() << std::endl;
//     std::cout << "Neighbors for T3: " << neighborsT3.size() << std::endl;
// }

// Function to calculate mean
double calculateMean(const std::vector<double>& data) {
    double sum = 0.0;
    for (double val : data) {
        sum += val;
    }
    return sum / data.size();
}

// Function to calculate standard deviation
double calculateStdDev(const std::vector<double>& data) {
    double mean = calculateMean(data);
    double variance = 0.0;
    for (double val : data) {
        variance += std::pow(val - mean, 2);
    }
    variance /= data.size();
    return std::sqrt(variance);
}

// Function to generate a textual histogram
void generateHistogram(const std::vector<double>& data, int bin_count) {
    double min_val = *std::min_element(data.begin(), data.end());
    double max_val = *std::max_element(data.begin(), data.end());
    double bin_width = (max_val - min_val) / bin_count;
    
    std::map<int, int> histogram;
    for (double val : data) {
        int bin = (val - min_val) / bin_width;
        histogram[bin]++;
    }
    
    std::cout << "Histogram:" << std::endl;
    for (const auto& bin : histogram) {
        std::cout << std::setw(4) << min_val + bin.first * bin_width << " - " << min_val + (bin.first + 1) * bin_width << " : ";
        for (int i = 0; i < bin.second; ++i) {
            std::cout << "*";
        }
        std::cout << std::endl;
    }
}

// Function to calculate mean and standard deviation
std::pair<double, double> meanAndStdDev(const std::vector<double>& vec) {
    double mean = std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
    double variance = std::accumulate(vec.begin(), vec.end(), 0.0, [mean](double acc, double val) {
        return acc + std::pow(val - mean, 2);
    }) / vec.size();
    double std_dev = std::sqrt(variance);
    return {mean, std_dev};
}

double variance(std::vector<double>& vec) {
    double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
    double mean = sum / vec.size();
    double sq_sum = std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
    return (sq_sum / vec.size()) - mean * mean;
}

// Calculate the variance of dot products of triangle normals
double calculateVarianceOfDotProducts(const std::vector<Triangle>& triangles) {
    std::vector<double> dotProducts;
    for (size_t i = 0; i < triangles.size(); ++i) {
        for (size_t j = i + 1; j < triangles.size(); ++j) {
            Vector3D normal_i = triangles[i].normal;
            Vector3D normal_j = triangles[j].normal;
            double dotProduct = normal_i.dot(normal_j);  // Assuming you have a dot function in Vector3D class
            dotProducts.push_back(dotProduct);
        }
    }

    return variance(dotProducts);
}

bool areNormalsEqual(const Vector3D& normal1, const Vector3D& normal2) {
    // You may need to adjust the threshold based on your data
    double threshold = 0.0001;
    return std::abs(normal1.x - normal2.x) < threshold &&
           std::abs(normal1.y - normal2.y) < threshold &&
           std::abs(normal1.z - normal2.z) < threshold;
}

bool areTrianglesCoplanar(const std::vector<Triangle>& triangles) {
    if (triangles.size() < 2) return true;  // Less than 2 triangles, considered coplanar

    // Access the pre-computed normal of the first triangle
    Vector3D firstNormal = triangles[0].normal;
    
    for (size_t i = 1; i < triangles.size(); ++i) {
        // Access the pre-computed normal of the current triangle
        Vector3D currentNormal = triangles[i].normal;
        
        if (!areNormalsEqual(firstNormal, currentNormal)) return false;
    }
    return true;
}

bool hasDuplicateVertices(const std::vector<Vector3D>& vertices) {
    std::unordered_set<Vector3D> vertexSet(vertices.begin(), vertices.end());
    return vertexSet.size() != vertices.size();
}

int main() {
    std::cout << "Starting program..." << std::endl;
    // std::vector<Triangle> allTriangles = {};  // Fill with triangles from your data source
    std::vector<Vector3D> allVertices;    
    double tau_S =.01;  // Threshold for shape similarity
    // std::string inputfile = "/home/kingpin/Documents/blender/church.obj";
    // std::string inputfile = "/home/kingpin/Documents/blender/burger.obj";
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
    std::cout << "Total triangles: " << allTriangles.size() << std::endl;
    std::cout << "Total vertices: " << allVertices.size() << std::endl;


    if (areTrianglesCoplanar(allTriangles)) {
        std::cerr << "Triangles are coplanar!" << std::endl;
    }

    if (hasDuplicateVertices(allVertices)) {
        std::cerr << "Found duplicate vertices!" << std::endl;
    }

    // exit(1);

    // Step 3: Data Analysis
    std::vector<double> areas;
    std::vector<double> edgeLengths;
    std::vector<double> minEdges;
    std::vector<double> maxEdges;
    std::vector<double> dotProducts;

    TriangleAdjacency triangleAdjacency;  // Populate this object with your adjacency data

    for (const Triangle& triangle : allTriangles) {
        triangleAdjacency.addTriangle(triangle);
    }

    for(const Triangle& t : allTriangles) {
        areas.push_back(t.area);
        minEdges.push_back(t.e_min);
        maxEdges.push_back(t.e_max);
        edgeLengths.push_back(t.e_max);

        auto adjacentTriangles = triangleAdjacency.getAdjacentTriangles(t);
        
        for (const Triangle& adjacentTriangle : adjacentTriangles) {
            Vector3D normal1 = t.normal;
            Vector3D normal2 = adjacentTriangle.normal;

            // Normalize the normals
            normal1.normalize();
            normal2.normalize();

            // Calculate the dot product
            double dotProduct = normal1.dot(normal2);

            dotProducts.push_back(dotProduct);
        }
    }

    // Calculate variances
    double varArea = variance(areas);
    double varMinEdge = variance(minEdges);
    double varMaxEdge = variance(maxEdges);
    double varDotProduct = variance(dotProducts);  // Calculate variance of dot products

    std::cout << "Variance of dot products: " << varDotProduct << std::endl;

    // Calculate weights
    // double sumVar = varArea + varMinEdge + varMaxEdge + varDotProduct;
    double sumVarClusters = varArea + varMinEdge + varMaxEdge;
    double wAClusters = varArea / sumVarClusters;
    double wLClusters = varMaxEdge / sumVarClusters;
    double wSClusters = varMinEdge / sumVarClusters;

    double sumVarComponents = varArea + varMinEdge + varMaxEdge + varDotProduct;
    double wAComponents = varArea / sumVarComponents;
    double wLComponents = varMaxEdge / sumVarComponents;
    double wSComponents = varMinEdge / sumVarComponents;
    double wNComponents = varDotProduct / sumVarComponents;  // Weight corresponding to dot product variance
    
    std::cout << "wA: " << wAClusters << ", wL: " << wLClusters << ", wS: " << wSClusters << std::endl;
    std::cout << "wAComponents: " << wAComponents << ", wLComponents: " << wLComponents << ", wSComponents: " << wSComponents << ", wNComponents: " << wNComponents << std::endl;
    std::vector<double> weightClusters = {wAClusters, wLClusters, wSClusters, 0};
    // std::vector<double> weightComponents = {wAComponents, wLComponents, wSComponents, wNComponents};
    std::vector<double> weightComponents = {wAComponents, wLComponents, wSComponents, wNComponents};
    // std::vector<double> weights = {wA, wL, wS, 0};

    // Step 2: Calculate Descriptive Statistics
    auto [meanArea, stdDevArea] = meanAndStdDev(areas);
    auto [meanEdge, stdDevEdge] = meanAndStdDev(edgeLengths);
    
    std::cout << "Mean Area: " << meanArea << ", Std Dev Area: " << stdDevArea << std::endl;
    std::cout << "Mean Edge Length: " << meanEdge << ", Std Dev Edge Length: " << stdDevEdge << std::endl;

    std::sort(edgeLengths.begin(), edgeLengths.end());

    // 1. Median Edge Length
    double medianEdgeLength = (edgeLengths[edgeLengths.size()/2] + edgeLengths[(edgeLengths.size() - 1)/2]) / 2.0;
    
    // Step 3: Choose an Optimal Cell Size
    double optimalCellSize = medianEdgeLength;  // This is a simple example; you can use more complex logic based on your needs
    
    std::cout << "Optimal Cell Size: " << optimalCellSize << std::endl;

    std::cout << "First 10 triangles:" << std::endl;
    for (int i = 0; i < 10 && i < allTriangles.size(); i++) {
        std::cout << "Triangle " << i << ": v1(" << allTriangles[i].vec1.x << "," << allTriangles[i].vec1.y << "," << allTriangles[i].vec1.z 
                << ") v2(" << allTriangles[i].vec2.x << "," << allTriangles[i].vec2.y << "," << allTriangles[i].vec2.z 
                << ") v3(" << allTriangles[i].vec3.x << "," << allTriangles[i].vec3.y << "," << allTriangles[i].vec3.z << ")" << std::endl;
    }

    // testSpatialHash(allTriangles);

    std::cout << "Number of triangles to be clustered: " << allTriangles.size() << std::endl;

    std::vector<Cluster> clusters = createSearchSpaces(allTriangles, weightClusters, weightComponents , tau_S, optimalCellSize); 
    std::cout << "Finished createSearchSpaces." << std::endl;

    std::cout << "Program completed successfully." << std::endl;

    return 0;
}