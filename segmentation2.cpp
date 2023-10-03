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
#include <queue> // for std::priority_queue
#include <utility> // for std::pair
#include <algorithm> // for std::all_of
#include <cstdlib> // for std::exit
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
// #include "/tinyobj_loader_opt.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Polyhedron_3<K> Polyhedron;
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

class Hull {
private:
    std::vector<Point_3> points;
    Polyhedron P;

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


    // Add a triangle to the collection of points
    void add(const Triangle& triangle) {
        points.push_back(Point_3(triangle.vec1.x, triangle.vec1.y, triangle.vec1.z));
        points.push_back(Point_3(triangle.vec2.x, triangle.vec2.y, triangle.vec2.z));
        points.push_back(Point_3(triangle.vec3.x, triangle.vec3.y, triangle.vec3.z));
    }

    void computeHull() {
        
        if (points.size() < 4) {
            return;
        }
        try {
            CGAL::convex_hull_3(points.begin(), points.end(), P);
        } catch (const std::exception& e) {
            std::cerr << "Exception: " << e.what() << std::endl;
        }

    }


    double volume() const {
        if (P.is_empty()) {
            // Log or handle the error. The polyhedron is empty.
            return 0.0;
        }

        double total_volume = 0.0;

        for (auto facet = P.facets_begin(); facet != P.facets_end(); ++facet) {
            auto h = facet->halfedge();
            
            Point_3 A = h->vertex()->point();
            Point_3 B = h->next()->vertex()->point();
            Point_3 C = h->opposite()->vertex()->point();

            // Compute Q_F, we can take A as a point on face F
            Point_3 Q_F = A;

            // Compute the normal of the face using the correct vector type
            K::Vector_3 N_F = CGAL::normal(A, B, C);

            // Compute area of triangle
            double areaF = std::sqrt(N_F.squared_length()) * 0.5;
            N_F = N_F / std::sqrt(N_F.squared_length()); // Normalize the normal

            total_volume += (Q_F - CGAL::ORIGIN) * N_F * areaF;
        }
        //ENABLE DEBUG
        // std::cout << "Computed volume: " << std::abs(total_volume) / 3.0 << std::endl;
        return std::abs(total_volume) / 3.0;
    }

};

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


    void updateHull() {
        // Only update the hull if we have enough data to create it.
        if (!triangles.empty()) {
            convexHull = Hull(triangles);
        }
    }

    bool isValid() const {
        // std::cout << "[Debug] Checking validity for component with triangles:" << std::endl;
        // for (const auto& triangle : triangles) {
        //     std::cout << triangle.toString() << std::endl;  // Assumes Triangle has a good ostream operator
        // }
        bool contiguous = areTrianglesContiguous(triangles);
        // std::cout << "[Debug] Is contiguous: " << std::boolalpha << contiguous << std::endl;
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


// Function to calculate the weight of a component type.
double getWeightOfComponent(const std::unordered_map<std::string, int>& typeCounts, const std::string& typeString) {
    int frequency = typeCounts.count(typeString) ? typeCounts.at(typeString) : 1; // default frequency is 1 to avoid division by zero
    return 1.0 / std::pow(frequency, 2); // This is |Zi|^2 in the denominator
}


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
    std::set<int> toRemove;  // To keep track of indices to remove

    // Debug: Count total triangles before merging
    int totalTrianglesBefore = 0;
    for (const auto& component : components) {
        totalTrianglesBefore += component.triangles.size();
    }
    std::cout << "Total triangles before merging: " << totalTrianglesBefore << std::endl;

    for (int i = 0; i < components.size(); ++i) {
        if (toRemove.count(i) > 0) {
            continue;  // Skip already removed components
        }

        Hull hull_i = Hull(components[i].triangles); 
        double volume_i = hull_i.volume();

        for (int j = i + 1; j < components.size(); ++j) {
            if (toRemove.count(j) > 0) {
                continue;  // Skip already removed components
            }

            Hull hull_j = Hull(components[j].triangles);  // Computing the hull only once for efficiency
            double volume_j = hull_j.volume();

            // Check for type and significant overlap
            if (components[i].type == components[j].type && components[i].overlapsSignificantly(components[j])) {
                if (attemptToMergeComponents(components[i], components[j], adjacency)) {
                    for (const auto& triangle : components[i].triangles) {
                        std::cout << triangle.toString() << std::endl;
                    }

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


double Stitj(const Triangle& ti, const Triangle& tj, const std::vector<double>& weights) {
    double wA = weights[0], wL = weights[1], wS = weights[2], wN = weights[3];
    
    double areaDenominator = std::max(ti.area, tj.area);
    double maxLengthDenominator = ti.e_max + tj.e_max;
    double minLengthDenominator = ti.e_min + tj.e_min;

    double areaDifference = (areaDenominator != 0.0) ? wA * std::abs(ti.area - tj.area) / areaDenominator : 0.0;
    double maxLengthDifference = (maxLengthDenominator != 0.0) ? wL * std::abs(ti.e_max - tj.e_max) / maxLengthDenominator : 0.0;
    double minLengthDifference = (minLengthDenominator != 0.0) ? wS * std::abs(ti.e_min - tj.e_min) / minLengthDenominator : 0.0;
    double normalDifference = wN * (1 - (ti.normal.x * tj.normal.x + ti.normal.y * tj.normal.y + ti.normal.z * tj.normal.z)) / 2.0;

    double totalSum = areaDifference + maxLengthDifference + minLengthDifference + normalDifference;

    // std::cout << "Area Percentage " << (areaDifference/totalSum) * 100 << "%" << std::endl;
    // std::cout << "Max Edge Percentage " << (maxLengthDifference/totalSum) * 100 << "%" << std::endl;
    // std::cout << "Min Edge Percentage " << (minLengthDifference/totalSum) * 100 << "%" << std::endl;

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

std::vector<Cluster> initialClusteringByShapeSimilarity(const std::vector<Triangle>& triangles, const std::vector<double>& weights, double tau_S) {
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


std::vector<std::vector<std::pair<Triangle, double>>> getCandidateTrianglesForEachComponent(
    const std::vector<Component>& tempComponents, 
    const std::unordered_set<Triangle>& remainingTriangles,  
    size_t topN = 3  
) {
    std::vector<std::vector<std::pair<Triangle, double>>> candidateTrianglesForComponents;
    std::unordered_set<Triangle> alreadyChosenTriangles;

    for (const Component& component : tempComponents) {
        std::priority_queue<std::pair<double, Triangle>> candidatePriorityQueue;
        
        // Using .volume() directly, since it's not an optional anymore.
        double originalVolume = component.convexHull.volume();

        for (const Triangle& triangle : remainingTriangles) {
            if (alreadyChosenTriangles.find(triangle) != alreadyChosenTriangles.end()) {
                continue;
            }

            Component tempComponent = component;
            tempComponent.addTriangle(triangle);  // This will update the hull too

            if (!tempComponent.isValid()) {
                continue;  // Skip non-contiguous triangles
            }

            // Using .volume() directly
            double newVolume = tempComponent.convexHull.volume();

            double volumeChange = std::abs(newVolume - originalVolume);

            // We use negative volumeChange because priority_queue sorts in descending order
            candidatePriorityQueue.push({-volumeChange, triangle});
        }

        std::vector<std::pair<Triangle, double>> topCandidates;

        for (size_t i = 0; i < topN && !candidatePriorityQueue.empty(); ++i) {
            auto scoreTrianglePair = candidatePriorityQueue.top();
            candidatePriorityQueue.pop();

            double score = -scoreTrianglePair.first;
            Triangle triangle = scoreTrianglePair.second;

            topCandidates.push_back({triangle, score});
            alreadyChosenTriangles.insert(triangle);
        }

        candidateTrianglesForComponents.push_back(topCandidates);
    }

    return candidateTrianglesForComponents;
}


double calculateConvexHullChange(const Component& component, const Triangle& candidateTriangle) {
    Hull originalHull(component.triangles);
    double originalVolume = originalHull.volume();

    // Create a new Hull object by first copying the existing component
    Component tempComponent = component;

    // Then add the new triangle to this temporary component
    tempComponent.addTriangle(candidateTriangle);  // Assuming you have a way to add a triangle to a Component

    // Compute the hull for this new component
    Hull newHull(tempComponent.triangles);
    newHull.computeHull();  // I'm assuming you'll need to recompute the hull here
    
    double newVolume = newHull.volume();
    return std::abs(newVolume - originalVolume);  // Returns absolute change in volume
}

double calculateGlobalRepetitionFactor(const std::vector<double>& globalVolumeChangeHistory, const std::vector<double>& minConvexHullChanges) {
    // If there's not enough data, return a neutral factor
    if (globalVolumeChangeHistory.empty() || minConvexHullChanges.empty()) return 1.0;

    // Calculate the mean and variance of the historical volume changes
    double mean = std::accumulate(globalVolumeChangeHistory.begin(), globalVolumeChangeHistory.end(), 0.0) / globalVolumeChangeHistory.size();
    double sq_sum = std::inner_product(globalVolumeChangeHistory.begin(), globalVolumeChangeHistory.end(), globalVolumeChangeHistory.begin(), 0.0);
    double variance = (sq_sum / globalVolumeChangeHistory.size()) - (mean * mean);

    // If variance is zero, the system is perfectly repetitive
    if (variance == 0) return 1.0;

    // Calculate the mean of current minimum changes
    double currentMean = std::accumulate(minConvexHullChanges.begin(), minConvexHullChanges.end(), 0.0) / minConvexHullChanges.size();
  
    // Calculate a repetition factor based on how close the current mean is to the historical mean, and how low the variance is
    double distanceToMean = std::abs(currentMean - mean);
    double repetitionFactor = std::exp(-distanceToMean / variance);

    return repetitionFactor;
}

std::vector<std::vector<Triangle>> synchronizeGrowth(
    const std::vector<std::vector<std::pair<double, Triangle>>>& adjustedCandidateTriangles,
    const std::vector<double>& minConvexHullChanges,
    const std::vector<double>& componentWeights,  // Updated this line to match the new type
    const std::vector<double>& globalVolumeChangeHistory,
    bool weightImpactAlreadyApplied = false, // New flag
    bool getTopN = true,
    size_t N = 1)
{
    if (adjustedCandidateTriangles.empty()) {
        return {};
    }

    std::vector<std::vector<Triangle>> bestTrianglesForComponents(adjustedCandidateTriangles.size());

    double globalMinChange = *std::min_element(minConvexHullChanges.begin(), minConvexHullChanges.end());
    
    // Calculate global repetition factor
    double globalRepetitionFactor = calculateGlobalRepetitionFactor(globalVolumeChangeHistory, minConvexHullChanges);

    for (size_t i = 0; i < adjustedCandidateTriangles.size(); ++i) {
        std::vector<Triangle> topTriangles;
        std::priority_queue<std::pair<double, Triangle>> pq;

        for (const auto& candidatePair : adjustedCandidateTriangles[i]) {
            double change = candidatePair.first;
            const Triangle& candidate = candidatePair.second;

            change *= globalRepetitionFactor;  // Apply global repetition factor

            pq.push({-std::abs(change - globalMinChange), candidate});
        }

        while (!pq.empty() && topTriangles.size() < N) {
            topTriangles.push_back(pq.top().second);
            pq.pop();
        }
        bestTrianglesForComponents[i] = topTriangles;
    }

    return bestTrianglesForComponents;
}



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


std::vector<Triangle> findLargestDisjointSimilarTriangles(const std::vector<Triangle>& triangles, 
                                                          std::vector<double>& weights, double tau_S) {
    std::vector<Triangle> seeds;

    if (triangles.empty()) {
        std::cerr << "Error: No triangles to find seeds from.\n";
        return seeds;
    }

    // Sort triangles by area
    std::vector<Triangle> sortedTriangles = triangles;
    std::sort(sortedTriangles.begin(), sortedTriangles.end(), [](const Triangle& a, const Triangle& b) {
        return a.area < b.area;
    });


    // Set weights, with wN not equal to 0
    double wN = .05;
    weights[3] = wN;

    for (const Triangle& triangle : sortedTriangles) {
        bool isSimilar = false;
        for (const Triangle& seed : seeds) {
            if (Stitj(triangle, seed, weights) <= tau_S) { // Close to 0, meaning triangles are similar
                isSimilar = true;
                break;
            }
        }
        if (!isSimilar) {
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
                                std::vector<double>& weights,
                                double tau_S) {
    return findLargestDisjointSimilarTriangles(
        std::vector<Triangle>(remainingTriangles.begin(), remainingTriangles.end()),
        weights,
        tau_S
    );                           
}

std::vector<Component> initializeTempComponents(const std::vector<Triangle>& seeds,
                                                std::unordered_set<Triangle>& remainingTriangles) {
    std::vector<Component> tempComponents;
    for (const auto& seed : seeds) {
        if (remainingTriangles.erase(seed) != 0) {
            Component newComponent;  // Should work now with the default constructor
            newComponent.addTriangle(seed);  // This will also update the Hull
            tempComponents.push_back(newComponent);
        }
    }
    return tempComponents;
}

double calculateWeight(int repetitionCount) {
    return 1.0 / (std::pow(repetitionCount, 2));
}


bool isValidWeight(double newWi, double oldWi) {
    std::cout << "New weight: " << newWi << std::endl;
    std::cout << "Old weight: " << oldWi << std::endl;

    // We only check for less than as per the paper's instruction to lower the weight
    return (newWi < oldWi);
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

enum StallReason {
    None,
    BestTriangleUsed,
    NoBestTriangles
};

bool applyGreedyAlgorithm(
    std::vector<Component>& tempComponents,
    std::unordered_set<Triangle>& remainingTriangles,
    std::unordered_set<Triangle>& usedTriangles,
    std::vector<int>& repetitionCounts,
    std::vector<double>& componentWeights,
    std::vector<double>& globalVolumeChangeHistory,
    std::vector<std::stack<Triangle>>& issueTrianglesStacks,
    bool successfulRun,
    int componentIndex
) {

    static std::unordered_set<std::string> triedPaths;
    // Generate serialized key based on the current state of tempComponents and componentIndex
    std::string currentPathId = serializePathVector(tempComponents);
    std::string serializedKey = currentPathId + "|" + std::to_string(componentIndex);

    // Check if this path has been tried before
    if (triedPaths.find(serializedKey) != triedPaths.end()) {
        // This path has already been tried; you may want to return false or handle this accordingly
        return false;
    }

    // If not, proceed with the existing logic and add this path to triedPaths
    triedPaths.insert(serializedKey);
    // Now apply the greedy algorithm again but avoid the issue triangles
    auto candidateTriangles = getCandidateTrianglesForEachComponent(tempComponents, remainingTriangles);  // Assuming this function exists

    std::vector<double> minConvexHullChanges(tempComponents.size(), std::numeric_limits<double>::infinity());
    std::vector<std::vector<std::pair<double, Triangle>>> adjustedCandidateTriangles(tempComponents.size());

    // Modified loop for considering both primary and alternate candidates
    for (size_t i = 0; i < tempComponents.size(); ++i) {
        int repetitionCount = repetitionCounts[i];
        double weightImpact = componentWeights[i];

        for (const auto& candidatePair : candidateTriangles[i]) {
            const Triangle& candidate = candidatePair.first;
            double change = calculateConvexHullChange(tempComponents[i], candidate);
            change *= weightImpact;
            minConvexHullChanges[i] = std::min(minConvexHullChanges[i], change);
            adjustedCandidateTriangles[i].emplace_back(change, candidate);
        }
    }

    std::cout << "MIN CONVEX HULL CHANGE: " << minConvexHullChanges.size() << std::endl;

    size_t N = 3;  // Number of top candidates to retrieve. 3 is arbitrary; set your desired number.
    auto bestTrianglesForComponents = synchronizeGrowth(adjustedCandidateTriangles, minConvexHullChanges, componentWeights, globalVolumeChangeHistory, true, N);
    std::cout << "Best triangles for components size: " << bestTrianglesForComponents.size() << std::endl;

    for (const Triangle& alternativeTriangle : bestTrianglesForComponents[componentIndex]) {
        // Skip the problematic triangles
        if (!issueTrianglesStacks[componentIndex].empty() && issueTrianglesStacks[componentIndex].top() == alternativeTriangle) {
            continue;
        }
        
        if (usedTriangles.find(alternativeTriangle) == usedTriangles.end()) {
                // Your existing validation logic here, including weight and convex hull checks
            int oldZi = repetitionCounts[componentIndex]; // Store the current repetition count
            int hypotheticalNewZi = oldZi + 1;
            double oldWi = calculateWeight(oldZi);
            double hypotheticalNewWi = calculateWeight(hypotheticalNewZi);

            if (isValidAfterAddingTriangle(tempComponents[componentIndex], alternativeTriangle) && 
                isValidWeight(hypotheticalNewWi, oldWi)) {
                    tempComponents[componentIndex].addTriangle(alternativeTriangle);  // Add the triangle to the component
                    repetitionCounts[componentIndex] = hypotheticalNewZi; // Update the repetition count to the hypothetical new value
                    componentWeights[componentIndex] = hypotheticalNewWi;  // Update the weight to the hypothetical new value

                    std::cout << "IS VALID AFTER ADDING TRIANGLE " << true << std::endl;

                    double convexHullChange = minConvexHullChanges[componentIndex];
                    std::cout << "Convex hull change for component " << componentIndex << ": " << convexHullChange << std::endl;
                    
                    remainingTriangles.erase(alternativeTriangle);
                    usedTriangles.insert(alternativeTriangle);
                    tempComponents[componentIndex].id = componentIndex;
                    globalVolumeChangeHistory.push_back(convexHullChange);
                    successfulRun = true;
            }
            else {
                successfulRun = false;
            }

        }
    }
    return successfulRun;
}

bool recursiveBackTracking(
    std::stack<std::vector<Component>>& globalBacktrackStack,
    std::stack<std::vector<Component>>& unhappyBacktrackStack,
    std::vector<std::stack<Triangle>>& issueTrianglesStacks,
    std::vector<StallReason> stallReasons,
    std::stack<int> componentIndexStack,
    std::vector<Component>& tempComponents,
    std::unordered_set<Triangle>& remainingTriangles,
    std::unordered_set<Triangle>& usedTriangles,
    std::vector<int>& repetitionCounts,
    std::vector<double>& componentWeights,
    std::vector<double>& globalVolumeChangeHistory,
    int depth)
{

    if(depth >= 3) {
        return false;
    }

    // Check if we need to backtrack to a good state
    if (!globalBacktrackStack.empty()) {
        tempComponents = globalBacktrackStack.top();  // Set tempComponents to the last good state
        globalBacktrackStack.pop();  // Remove the last good state from the stack

        // Optionally, clear unhappyBacktrackStack, issueTrianglesStacks, and other data structures here
    } else {
        return false;  // No good states to revert to
    }

    if (componentIndexStack.size() > 1) {
        // Full Backtrack
        //Maybe instead of popping the entire top stack, we can just pop the top element
        // New Logic: Identify components that have extra triangles (fulfilled but overgrew)
        std::vector<int> overgrewComponentIndices;
        // Logic to fill overgrewComponentIndices; this might involve comparing globalBacktrackStack and unhappyBacktrackStack.

        // New Logic: Remove the last added triangle for each overgrew component in tempComponents
        for (const int i : overgrewComponentIndices) {
            Triangle lastAddedTriangle = tempComponents[i].getLastTriangle();
            tempComponents[i].removeTriangle(lastAddedTriangle);

            // Update remainingTriangles and usedTriangles
            remainingTriangles.insert(lastAddedTriangle);
            usedTriangles.erase(lastAddedTriangle);
        }

        bool successfulRun = false; // Keep track if the run is successful

        // Remove the last added triangle from each specific component in tempComponents
        for (const int i : specificComponentIndices) {
            Triangle lastAddedTriangle = tempComponents[i].getLastTriangle();  // Assuming getLastTriangle() exists
            tempComponents[i].removeTriangle(lastAddedTriangle); // Remove from component
            remainingTriangles.insert(lastAddedTriangle);  // Put back into remainingTriangles
            usedTriangles.erase(lastAddedTriangle);  // Remove from usedTriangles

            successfulRun = applyGreedyAlgorithm(tempComponents, remainingTriangles, usedTriangles, 
            repetitionCounts, componentWeights, 
            globalVolumeChangeHistory, issueTrianglesStacks,
            successfulRun, i);
        }

        if (!successfulRun) {
            // Recursively call the backtracking function with incremented depth
            return recursiveBackTracking(
                globalBacktrackStack,
                unhappyBacktrackStack,
                issueTrianglesStacks,
                stallReasons,
                componentIndexStack,
                tempComponents,
                remainingTriangles,
                usedTriangles,
                repetitionCounts,
                componentWeights,
                globalVolumeChangeHistory,
                depth + 1  // Increment depth
            );
        }

        // If the run is successful, clear the issue triangles as they are no longer relevant
        if (successfulRun) {
            for (auto& stack : issueTrianglesStacks) {
                while (!stack.empty()) {
                    stack.pop();
                }
            }
            //remove the componentIndexStack so it doesn't go in the below while loop
            while (!componentIndexStack.empty()) {
                int componentIndex = componentIndexStack.top();
                componentIndexStack.pop();
            }
        } else {
            return false;
        }
    }


    if (componentIndexStack.size() == 1) {
        // Initialize variables for this iteration
        bool anyComponentGrown = false;
        bool allComponentsStalled = true;
        std::vector<Component> grownTempComponents;
        std::vector<Component> badlyGrownTempComponents;
        std::vector<bool> componentsStalled(tempComponents.size(), false);

        int failedComponentIndex = -1;
        if (!componentIndexStack.empty()) {
            failedComponentIndex = componentIndexStack.top();
            componentIndexStack.pop();  // Remove this index as it is now being processed
        }

        auto candidateTriangles = getCandidateTrianglesForEachComponent(tempComponents, remainingTriangles);

        std::cout << "Candidate triangles size: " << candidateTriangles.size() << std::endl;

        std::vector<double> minConvexHullChanges(tempComponents.size(), std::numeric_limits<double>::infinity());
        std::vector<std::vector<std::pair<double, Triangle>>> adjustedCandidateTriangles(tempComponents.size());

        // Modified loop for considering both primary and alternate candidates
        for (size_t i = 0; i < tempComponents.size(); ++i) {
            int repetitionCount = repetitionCounts[i];
            double weightImpact = componentWeights[i];

            for (const auto& candidatePair : candidateTriangles[i]) {
                const Triangle& candidate = candidatePair.first;
                double change = calculateConvexHullChange(tempComponents[i], candidate);
                change *= weightImpact;
                minConvexHullChanges[i] = std::min(minConvexHullChanges[i], change);
                adjustedCandidateTriangles[i].emplace_back(change, candidate);
            }
        }

        std::cout << "MIN CONVEX HULL CHANGE: " << minConvexHullChanges.size() << std::endl;

        size_t N = 3;  // Number of top candidates to retrieve. 3 is arbitrary; set your desired number.
        auto bestTrianglesForComponents = synchronizeGrowth(adjustedCandidateTriangles, minConvexHullChanges, componentWeights, globalVolumeChangeHistory, true, N);
        std::cout << "Best triangles for components size: " << bestTrianglesForComponents.size() << std::endl;
        for (size_t i = (failedComponentIndex == -1 ? 0 : failedComponentIndex); i < tempComponents.size(); ++i) {
            StallReason currentStallReason = stallReasons[i];

            if (currentStallReason == StallReason::BestTriangleUsed) {
                Triangle issueTriangle;
                if (!issueTrianglesStacks[i].empty()) {
                    issueTriangle = issueTrianglesStacks[i].top();  // Assuming the stack's top contains the latest issue triangle
                }

                for (const Triangle& alternativeTriangle : bestTrianglesForComponents[i]) {
                    // Skip the problematic triangle
                    if (issueTriangle == alternativeTriangle) {
                        continue;
                    }
                    // if (alreadyTriedTrianglesForComponents[i].find(alternativeTriangle) != alreadyTriedTrianglesForComponents[i].end()) {
                    //     continue;
                    // }

                    // Check if the alternative triangle is unused
                    if (usedTriangles.find(alternativeTriangle) == usedTriangles.end()) {
                         // Your existing validation logic here, including weight and convex hull checks
                        int oldZi = repetitionCounts[i]; // Store the current repetition count
                        int hypotheticalNewZi = oldZi + 1;
                        double oldWi = calculateWeight(oldZi);
                        double hypotheticalNewWi = calculateWeight(hypotheticalNewZi);

                        if (isValidAfterAddingTriangle(tempComponents[i], alternativeTriangle) && 
                            isValidWeight(hypotheticalNewWi, oldWi)) {
                                tempComponents[i].addTriangle(alternativeTriangle);  // Add the triangle to the component
                                repetitionCounts[i] = hypotheticalNewZi; // Update the repetition count to the hypothetical new value
                                componentWeights[i] = hypotheticalNewWi;  // Update the weight to the hypothetical new value

                                std::cout << "IS VALID AFTER ADDING TRIANGLE " << true << std::endl;

                                double convexHullChange = minConvexHullChanges[i];
                                std::cout << "Convex hull change for component " << i << ": " << convexHullChange << std::endl;
                                
                                remainingTriangles.erase(alternativeTriangle);
                                anyComponentGrown = true;
                                usedTriangles.insert(alternativeTriangle);
                                // componentIndexStack.push(i);
                                tempComponents[i].id = i;
                                grownTempComponents.push_back(tempComponents[i]);
                                globalVolumeChangeHistory.push_back(convexHullChange);
                            }
                    }
                }
            }
            else if (currentStallReason == StallReason::NoBestTriangles) {
                // Case 2: No Best Triangles
                // Backtrack to a state where the component had fewer triangles.
                // Remove the last added triangle.
                // Reapply the greedy algorithm.
                // TODO: Implement your logic here
            }

            // ... (Your existing logic for selecting and adding triangles)

        }

        bool anyComponentStalled = std::any_of(componentsStalled.begin(), componentsStalled.end(), [](bool v) { return v; });

        if (anyComponentStalled) {
            globalBacktrackStack.push(grownTempComponents);
            unhappyBacktrackStack.push(badlyGrownTempComponents);

            // Recursive call
            if (recursiveBackTracking(
                    globalBacktrackStack, unhappyBacktrackStack, issueTrianglesStacks, stallReasons,
                    componentIndexStack, tempComponents, remainingTriangles,
                    usedTriangles,
                    repetitionCounts, componentWeights, globalVolumeChangeHistory, depth + 1)) {
                return true;
            }

            // If we reached this point, recursive backtracking failed
            // Pop the last state from the global backtrack stack
            if (!globalBacktrackStack.empty()) {
                tempComponents = globalBacktrackStack.top();
                globalBacktrackStack.pop();
            }
        }
    }
    return false;
}


bool recursiveBackTracking(
    std::stack<std::vector<Component>>& globalBacktrackStack,
    std::stack<std::vector<Component>>& unhappyBacktrackStack,
    std::vector<std::stack<Triangle>>& issueTrianglesStacks,
    std::vector<StallReason> stallReasons,
    std::stack<int> componentIndexStack,
    std::vector<Component>& tempComponents,
    std::unordered_set<Triangle>& remainingTriangles,
    std::unordered_set<Triangle>& usedTriangles,
    std::vector<int>& repetitionCounts,
    std::vector<double>& componentWeights,
    std::vector<double>& globalVolumeChangeHistory,
    int depth)
{
    std::cout << "Entering recursiveBackTracking with depth: " << depth << std::endl;

    // Right before returning false for max depth
    if (depth >= 3) {
        std::cout << "Maximum depth reached at depth: " << depth << ". Backtracking failed." << std::endl;
        return false;
    }

    // Right before each recursive call
    std::cout << "Making recursive call with incremented depth: " << (depth + 1) << std::endl;

    // Right before returning false for each failure condition
    std::cout << "Returning false due to condition XYZ at depth: " << depth << std::endl;

    std::cout << "Before serialization, vector size: " << tempComponents.size() << std::endl;

    std::string currentPathId = serializePathVector(tempComponents);
    // std::cout << "Current path ID: " << currentPathId << std::endl;

    std::cout << "BEFORE POPPING: "<< std::endl;
    
    inspectStack(globalBacktrackStack);
    std::cout << "Unhappy backtrack stack: " << std::endl;
    inspectStack(unhappyBacktrackStack);
    for (size_t i = 0; i < issueTrianglesStacks.size(); ++i) {
        std::cout << "Issue triangles for component " << i << ": " << std::endl;
        
        // Note that printStack will empty the stack,
        // so make a copy if you need to keep the original intact
        std::stack<Triangle> tempStack = issueTrianglesStacks[i];
        printStack(tempStack);
    }
    std::cout << "--------Component index stack:-----------" << std::endl;
    printIndexStack(componentIndexStack);
    std::vector<Component> tempComponent;
    if (!globalBacktrackStack.empty()) {
        tempComponent = globalBacktrackStack.top();
        globalBacktrackStack.pop();
    }
    // std::cout << "AFTER POPPING: "<< std::endl;
    // inspectStack(globalBacktrackStack);
    std::cout << "Temp components size: " << tempComponent.size() << std::endl;
    for(const auto& comp : tempComponent) {
        std::cout << "Components AREEEEEEEEEe: " << comp.toString() << std::endl;
    }


    int multipleStalls = std::count_if(issueTrianglesStacks.begin(), issueTrianglesStacks.end(), [](const std::stack<Triangle>& s) {
        return !s.empty();
    });

    std::cout << multipleStalls << " components stalled." << std::endl;
    

    std::cout << "Temp components size: " << tempComponents.size() << std::endl;
    
    auto candidateTriangles = getCandidateTrianglesForEachComponent(tempComponents, remainingTriangles);

    std::cout << "Candidate triangles size: " << candidateTriangles.size() << std::endl;

    std::vector<double> minConvexHullChanges(tempComponents.size(), std::numeric_limits<double>::infinity());
    std::vector<std::vector<std::pair<double, Triangle>>> adjustedCandidateTriangles(tempComponents.size());

    // Modified loop for considering both primary and alternate candidates
    for (size_t i = 0; i < tempComponents.size(); ++i) {
        int repetitionCount = repetitionCounts[i];
        double weightImpact = componentWeights[i];

        for (const auto& candidatePair : candidateTriangles[i]) {
            const Triangle& candidate = candidatePair.first;
            double change = calculateConvexHullChange(tempComponents[i], candidate);
            change *= weightImpact;
            minConvexHullChanges[i] = std::min(minConvexHullChanges[i], change);
            adjustedCandidateTriangles[i].emplace_back(change, candidate);
        }
    }

    std::cout << "MIN CONVEX HULL CHANGE: " << minConvexHullChanges.size() << std::endl;

    size_t N = 3;  // Number of top candidates to retrieve. 3 is arbitrary; set your desired number.
    auto bestTrianglesForComponents = synchronizeGrowth(adjustedCandidateTriangles, minConvexHullChanges, componentWeights, globalVolumeChangeHistory, true, N);
    std::cout << "Best triangles for components size: " << bestTrianglesForComponents.size() << std::endl;
    for (size_t i = 0; i < tempComponents.size(); ++i) {
        for (const Triangle& triangle : bestTrianglesForComponents[i]) {

            std::vector<Component> oldComponents = tempComponents;
            std::vector<int> oldRepetitionCounts = repetitionCounts;
            std::vector<double> oldComponentWeights = componentWeights;

            bool allValid = true;

            std::vector<double> localConvexHullChanges;

            for (size_t j = 0; j < tempComponents.size(); ++j) {
                int oldZi = repetitionCounts[j];
                int hypotheticalNewZi = oldZi + 1;

                double oldWi = componentWeights[j];
                double hypotheticalNewWi = calculateWeight(hypotheticalNewZi);

                double hypotheticalConvexHullChange = calculateConvexHullChange(tempComponents[j], bestTrianglesForComponents[j][0]);

                if (!isValidAfterAddingTriangle(tempComponents[j], bestTrianglesForComponents[j][0]) || 
                    !isValidWeight(hypotheticalNewWi, oldWi)) {
                    allValid = false;
                    break;
                }

                localConvexHullChanges.push_back(hypotheticalConvexHullChange);
            }

            if (allValid) {
                for (size_t j = 0; j < tempComponents.size(); ++j) {
                    tempComponents[j].addTriangle(bestTrianglesForComponents[j][0]);
                    repetitionCounts[j] = oldRepetitionCounts[j] + 1;
                    componentWeights[j] = calculateWeight(repetitionCounts[j]);
                    globalVolumeChangeHistory.push_back(minConvexHullChanges[j]); // Update the global volume change history
                }
                
                globalBacktrackStack.push(tempComponents);

                if (recursiveBackTracking(globalBacktrackStack, unhappyBacktrackStack, issueTrianglesStacks, stallReasons,
                                        componentIndexStack, tempComponents, remainingTriangles,
                                        usedTriangles,
                                        repetitionCounts, componentWeights,
                                        globalVolumeChangeHistory, depth + 1)) {
                    return true;
                }
            }

            // Revert the changes
            tempComponents = oldComponents;
            repetitionCounts = oldRepetitionCounts;
            componentWeights = oldComponentWeights;

            // Update the tried triangles and paths
            std::string newPathId = serializePathVector(tempComponents);
        }
    }

    return false;
}


//Does the seed grow with the temp component? like is the seed attached?
//Why are triangles shrinking so rapidly?
bool growComponentsSynchronously(
    std::vector<Component>& tempComponents,
    std::unordered_set<Triangle>& remainingTriangles,
    TriangleAdjacency& adjacency,
    std::unordered_map<std::string, int>& typeCounts,
    std::unordered_set<Triangle>& usedTriangles) {

    std::cout << "Initial number of components: " << tempComponents.size() << std::endl;

    // Existing variable declarations
    std::vector<int> repetitionCounts(tempComponents.size(), 0);  // Count repetitions for each component
    std::vector<double> componentWeights(tempComponents.size(), 1.0);  // Initial weights for each component
    std::vector<std::stack<std::vector<Triangle>>> backtrackCandidateStacks(tempComponents.size());
    std::vector<bool> componentsStalled(tempComponents.size(), false);
    bool anyComponentFailed = false;
    std::vector<double> globalVolumeChangeHistory;
    std::vector<std::vector<std::pair<double, Triangle>>> adjustedCandidateTriangles(tempComponents.size());
    std::stack<std::vector<Component>> globalBacktrackStack;  
    std::stack<std::vector<Component>> unhappyBacktrackStack;  
    std::vector<std::stack<Triangle>> issueTrianglesStacks(tempComponents.size());
    std::stack<int> componentIndexStack;
    std::vector<StallReason> stallReasons(tempComponents.size(), StallReason::None);



    for (size_t i = 0; i < tempComponents.size(); ++i) {
        repetitionCounts[i] = 1;
        componentWeights[i] = 1.0;
    }

    if (tempComponents.empty() || remainingTriangles.empty()) {
        return false;
    }
    bool anyComponentGrown = false;
    bool allComponentsStalled = false;
    bool backtrackingSuccessful = false;
    std::vector<Component> grownTempComponents;
    std::vector<Component> badlyGrownTempComponents;

    std::vector<Component> updatedComponents;  // Step 1
    

    while (true) {
        bool allComponentsSuccessfullyGrown = true;

        std::cout << "----------While loop reset for next iteration-----------" << std::endl;

        // New Code: Reset flags for each iteration
        anyComponentGrown = false;
        allComponentsStalled = false; // New variable
        
        std::fill(componentsStalled.begin(), componentsStalled.end(), false);

        // Retrieve primary candidate triangles
        auto candidateTriangles = getCandidateTrianglesForEachComponent(tempComponents, remainingTriangles);

        std::vector<double> minConvexHullChanges(tempComponents.size(), std::numeric_limits<double>::infinity());
        // Modified loop for considering both primary and alternate candidates
        for(size_t i = 0; i < tempComponents.size(); ++i) {
            int repetitionCount = repetitionCounts[i];
            double weightImpact = componentWeights[i];
            std::cout << "Component " << i << " repetition count : " << repetitionCount << std::endl;
            std::cout << "Component " << i << " weight impact: " << weightImpact << std::endl;

            for(const auto& candidatePair : candidateTriangles[i]) {
                const Triangle& candidate = candidatePair.first;
                double change = calculateConvexHullChange(tempComponents[i], candidate);
                change *= weightImpact;
                minConvexHullChanges[i] = std::min(minConvexHullChanges[i], change);
                
                adjustedCandidateTriangles[i].emplace_back(change, candidate);
            }
        }
        
        size_t N = 3;  // Number of top candidates to retrieve. 3 is arbitrary; set your desired number.
        auto bestTrianglesForComponents = synchronizeGrowth(adjustedCandidateTriangles, minConvexHullChanges, componentWeights, globalVolumeChangeHistory, true, N);
        std::cout << "Best triangles for components size: " << bestTrianglesForComponents.size() << std::endl;
        for (size_t i = 0; i < tempComponents.size(); ++i) {
            std::vector<Triangle> viableTriangles = bestTrianglesForComponents[i];
            std::cout << "Viable triangles for component " << i << ": " << viableTriangles.size() << std::endl;
            if (!viableTriangles.empty()) {
                backtrackCandidateStacks[i].push(viableTriangles);
            }
        }

        for (size_t i = 0; i < tempComponents.size(); ++i) {

            if (bestTrianglesForComponents[i].empty()) {
                std::cout << "No best triangles for component " << i << std::endl;
                stallReasons[i] = StallReason::NoBestTriangles;
                allComponentsSuccessfullyGrown = false;
                tempComponents[i].id = i;
                badlyGrownTempComponents.push_back(tempComponents[i]);
                componentIndexStack.push(i);
                componentsStalled[i] = true;
                continue;
            }

            Triangle bestTriangle = bestTrianglesForComponents[i][0];
            std::cout << "Best triangle for component " << i << ": " << bestTriangle.toString() << std::endl;
            
            if (usedTriangles.find(bestTriangle) != usedTriangles.end()) {
                std::cout << "Best triangle for component " << i << " already used." << std::endl;
                std::cout << "Best triangle: " << bestTriangle.toString() << std::endl;
                stallReasons[i] = StallReason::BestTriangleUsed;
                issueTrianglesStacks[i].push(bestTriangle);
                allComponentsSuccessfullyGrown = false;
                componentsStalled[i] = true;
                tempComponents[i].id = i;
                badlyGrownTempComponents.push_back(tempComponents[i]);
                componentIndexStack.push(i);
                continue;
            }

            int oldZi = repetitionCounts[i]; // Store the current repetition count
            int hypotheticalNewZi = oldZi + 1; // This would be the new repetition count if we add the triangle

            double oldWi = calculateWeight(oldZi);  // Calculate the old weight
            double hypotheticalNewWi = calculateWeight(hypotheticalNewZi);  // Pre-calculate the potential new weight
            if (isValidAfterAddingTriangle(tempComponents[i], bestTriangle) && 
                isValidWeight(hypotheticalNewWi, oldWi)) {  // Use the hypothetical new weight
                tempComponents[i].addTriangle(bestTriangle);  // Add the triangle to the component
                repetitionCounts[i] = hypotheticalNewZi; // Update the repetition count to the hypothetical new value
                componentWeights[i] = hypotheticalNewWi;  // Update the weight to the hypothetical new value

                std::cout << "IS VALID AFTER ADDING TRIANGLE " << std::endl;
                double convexHullChange = minConvexHullChanges[i];
                
                remainingTriangles.erase(bestTriangle);
                anyComponentGrown = true;
                usedTriangles.insert(bestTriangle);
                tempComponents[i].id = i;
                grownTempComponents.push_back(tempComponents[i]);
                globalVolumeChangeHistory.push_back(convexHullChange);

                updatedComponents.push_back(tempComponents[i]);  // Step 2
 
            } else {
                std::cout << "----------INVALID TRIANGLE----------" << std::endl;
                std::cout << "Component " << i << " failed to grow." << std::endl;
                tempComponents[i].id = i;
                badlyGrownTempComponents.push_back(tempComponents[i]);
                componentIndexStack.push(i);
                componentsStalled[i] = true;
                allComponentsStalled = true;
            }

            // Output debugging info
            std::cout << "Component " << i << " repetition count : " << hypotheticalNewZi << std::endl;
            std::cout << "Component " << i << " weight impact: " << hypotheticalNewWi << std::endl;
            std::cout << "newWi: " << hypotheticalNewWi << std::endl;
            std::cout << "oldWi: " << componentWeights[i] << std::endl;
            std::cout << "IS VALID AFTER ADDING TRIANGLE " << isValidAfterAddingTriangle(tempComponents[i], bestTriangle) << std::endl;
            std::cout << "IS VALID WEIGHT " << isValidWeight(hypotheticalNewWi, componentWeights[i]) << std::endl;
            std::cout << "Temp component " << i << " triangles after trying to grow: " << tempComponents[i].triangles.size() << std::endl;
        }
        if (allComponentsSuccessfullyGrown) {
            globalBacktrackStack.push(tempComponents);
        }
        tempComponents.clear();  // Step 3
        tempComponents = updatedComponents;
        updatedComponents.clear();  // Clear it for the next iteration

        for (size_t i = 0; i < tempComponents.size(); ++i) {
            if (componentsStalled[i]) {
                std::cout << "Component " << i << " is stalled." << std::endl;
                std::cout << "Backtrack candidates for component " << i << ": " << backtrackCandidateStacks[i].size() << std::endl;
            }
        }
        
        bool anyComponentStalled = std::any_of(componentsStalled.begin(), componentsStalled.end(), [](bool v) { return v; });

        if (anyComponentStalled) {
            unhappyBacktrackStack.push(badlyGrownTempComponents);
            std::cout << "-----BACKTRACKING NOW-----" << std::endl;
                bool result = recursiveBackTracking(
                    globalBacktrackStack, unhappyBacktrackStack, issueTrianglesStacks, stallReasons,
                    componentIndexStack, tempComponents, remainingTriangles, 
                    usedTriangles, repetitionCounts, componentWeights, globalVolumeChangeHistory, 0
                );

            if (!backtrackingSuccessful) {
                std::cout << "Maximum depth reached. Backtracking failed." << std::endl;
                anyComponentGrown = false;
            }
        }

        if (!anyComponentGrown) {
            break;
        }
    }
    return anyComponentGrown;
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


// void reseedComponents(
//     std::vector<Component>& tempComponents,
//     std::unordered_set<Triangle>& remainingTriangles,
//     std::vector<double>& weights,
//     double tau_S,
//     std::unordered_map<std::string, int>& typeCounts
// ) {
//     std::vector<Triangle> newSeeds = findSeeds(remainingTriangles, weights, tau_S);

//     // Update typeCounts based on the new seeds
//     for (const auto& seed : newSeeds) {
//         typeCounts[seed.type]++;  // Increment the count for this type
//     }

//     // You may also decide to sort newSeeds based on typeCounts, if desired
//     std::sort(newSeeds.begin(), newSeeds.end(), [&](const Triangle& a, const Triangle& b) {
//         return typeCounts[a.type] > typeCounts[b.type];
//     });

//     // Restore skippedTriangles back to remainingTriangles
//     for(const auto& skippedTriangle : skippedTriangles) {
//         remainingTriangles.insert(skippedTriangle);
//     }
//     skippedTriangles.clear();

//     // Re-initialize tempComponents with the new seeds
//     if(!newSeeds.empty()) {
//         tempComponents = initializeTempComponents(newSeeds, remainingTriangles);
//     }
// }

//Create a final components that the tempComponents from the first run is added to,
//I think reseed, restructures the tempComponents to start from scratch
//I dont know why: 
// Temp component triangles after growing: 5
// Temp component triangles after growing: 3
// Temp component triangles after growing: 6
// Temp component triangles after growing: 5
//Maybe they're not similar triangles that are growing? Possibly have to save tempComponents pre and post growth to obj file
//They're growing so unevenly. I would think this would be caught in the validation step but it's not
//I think it's because the triangles are being added to the tempComponents in a random order - maybe? (copilot suggestion)

std::vector<Component> randomizedGrowthOptimization(const Cluster& cluster,
                                                    TriangleAdjacency& adjacency,
                                                    std::vector<double>& weights,
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


    std::cout << "-----iteration: " << iteration << " starting----- " << std::endl;
    std::cout << "Number of triangles in cluster: " << cluster.triangles.size() << std::endl;

    while (!remainingTriangles.empty()) {
        std::vector<Component> currentIterationComponents;
        std::vector<Triangle> seeds = findSeeds(remainingTriangles, weights, tau_S);

        if (seeds.empty()) {
            std::cout << "No seeds found. Exiting." << std::endl;
            break;
        }

        std::vector<Component> tempComponents = initializeTempComponents(seeds, remainingTriangles);
        bool canGrow = true;

        while (canGrow) {
            std::cout << "Triangle count before growing in remainingTriangles: " << remainingTriangles.size() << std::endl;
            std::cout << "Triangle count in seeds: " << seeds.size() << std::endl;
            for(const auto& component : tempComponents) {
                std::cout << "Temp component triangles before growing: " << component.triangles.size() << std::endl;
            }
            canGrow = growComponentsSynchronously(
                tempComponents, remainingTriangles, adjacency,
                typeCounts, usedTriangles
            );

            std::cout << "--------DONE WITH GROW COMPONENTS SYNCHRONOUSLY---------" << std::endl;

            std::cout << "Triangle count after growing in remainingTriangles: " << remainingTriangles.size() << std::endl;
            std::cout << "Triangle count in seeds after growing: " << seeds.size() << std::endl;

            // Check contiguity of each component
            for (const auto& component : tempComponents) {
                std::cout << "Temp component triangles after growing: " << component.triangles.size() << std::endl;
                if (!areTrianglesContiguous(component.triangles)) {
                    std::cout << "Component became non-contiguous during growth." << std::endl;
                }
            }
            // Re-seed the remaining triangles here if necessary
            std::vector<Triangle> newSeeds = findSeeds(remainingTriangles, weights, tau_S);

            // if (!newSeeds.empty()) {
            //     reseedComponents(tempComponents, remainingTriangles, weights, tau_S, skippedTriangles, typeCounts);
            // }
        }

        // Adding components to finalComponents, but with uniqueness check
        for (const auto& comp : currentIterationComponents) {
            bool alreadyInFinal = std::any_of(finalComponents.begin(), finalComponents.end(),
                                            [&comp](const Component& existingComponent) { return existingComponent == comp; });
            if (!alreadyInFinal) {
                finalComponents.push_back(comp);
            }
        }

        std::cout << "Iteration " << iteration << " completed." << std::endl;
        iteration++;
    }
    // Before returning finalComponents, perform the merging operation
    std::cout << "Before merging: " << finalComponents.size() << std::endl;
    mergeOverlappingComponents(finalComponents, adjacency);
    std::cout << "Triangle count at the start: " << cluster.triangles.size() << std::endl;
    std::cout << "After merging: " << finalComponents.size() << std::endl;

    // Final validation
    if (!validateFinalState(finalComponents, cluster)) {
        std::cout << "Validation failed." << std::endl;
        exit(1);  // Exits the program with a status of 1
    }

    return finalComponents;
}

void optimizeAndCategorizeComponents(std::vector<Cluster>& allClusters,
                                     TriangleAdjacency& adjacency,
                                     std::vector<double>& weights,
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

        // Now update C and Z based on finalComponentsForCluster
        for (const Component& component : finalComponentsForCluster) {
            ComponentType type = identifyComponentType(component);
            if (C.find(type) != C.end()) {
                Z[type].push_back(component);
                componentWeights[type] = 1.0 / std::pow(Z[type].size(), 2);
            } else {
                C.insert(type);
                Z[type] = {component};
                componentWeights[type] = 1.0;
            }
        }
    }
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


std::vector<Cluster> createSearchSpaces(const std::vector<Triangle>& triangles, const std::vector<double>& weights, 
                                        double tau_S, double spatialHashCellSize) {    // 1. Initial clustering based on shape similarity
    std::cout << "1. Initial clustering based on shape similarity" << std::endl;
    std::vector<double> within_cluster_stitj;
    std::vector<double> between_cluster_stitj;

    // Your clustering function call
    // std::vector<Cluster> clusters = initialClusteringByShapeSimilarity(triangles, weights, tau_S, maxArea, maxEdgeMax, maxEdgeMin, within_cluster_stitj, between_cluster_stitj);

    std::vector<Cluster> clusters = initialClusteringByShapeSimilarity(triangles, weights, tau_S);

    std::cout << "Number of initial clusters: " << clusters.size() << std::endl;
    std::cout << "Initial clustering done. Number of clusters: " << clusters.size() << std::endl;
    
    // 3. Iterative merging and splitting based on adjacency similarity
    std::cout << "3. Iterative merging and splitting based on adjacency similarity" << std::endl;
    double tau_N = 0.5;  // Threshold for adjacency similarity
    refineClusters(clusters, tau_N, spatialHashCellSize);

    std::cout << "Number of clusters after refining: " << clusters.size() << std::endl;

    std::vector<Component> allFinalComponents;
    
    std::vector<double> weightsPostClustering = {.25, .25, .25, 0};

    // // 4. Randomized Growth Optimization
    std::cout << "4. Randomized Growth Optimization" << std::endl;
    TriangleAdjacency adjacencyInstance; // Assume you've created and possibly initialized this instance
    for (const Cluster& cluster : clusters) {
        for (const Triangle& triangle : cluster.triangles) {
            adjacencyInstance.addTriangle(triangle);
        }
    }

    optimizeAndCategorizeComponents(clusters, adjacencyInstance, weightsPostClustering, tau_S, allFinalComponents);
    
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

    // Step 3: Data Analysis
    std::vector<double> areas;
    std::vector<double> edgeLengths;
    std::vector<double> minEdges;
    std::vector<double> maxEdges;

    // std::vector<double> weights = {wA, wL, wS, 0};
    for(const Triangle& t : allTriangles) {
        areas.push_back(t.area);
        minEdges.push_back(t.e_min);  // Assuming e_min is the minimum edge length for the triangle
        maxEdges.push_back(t.e_max);  // Assuming e_max is the maximum edge length for the triangle
        edgeLengths.push_back(t.e_max);
    }

    // Calculate variances
    double varArea = variance(areas);
    double varMinEdge = variance(minEdges);
    double varMaxEdge = variance(maxEdges);

    // Calculate weights
    double sumVar = varArea + varMinEdge + varMaxEdge;
    double wA = varArea / sumVar;
    double wL = varMaxEdge / sumVar;
    double wS = varMinEdge / sumVar;
    
    std::cout << "wA: " << wA << ", wL: " << wL << ", wS: " << wS << std::endl;

    std::vector<double> weights = {wA, wL, wS, 0};
    // std::vector<double> weights = {.3, .3, .3, 0};

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

    std::vector<Cluster> clusters = createSearchSpaces(allTriangles, weights, tau_S, optimalCellSize); 
    std::cout << "Finished createSearchSpaces." << std::endl;

    // std::cout << "Total triangles using indices: " << allTriangleIndices.size() << std::endl;

    // std::string outputFilename = "simple_output.obj";
    // simpleWriteOBJ(allVertices, allTriangleIndices, outputFilename);

    // // 1. Create search spaces by clustering similar triangles.
    // std::cout << "Starting createSearchSpaces..." << std::endl;
    // allTriangles = indicesToTriangles(allTriangleIndices, allVertices);

    // std::cout << "First 10 vertices:" << std::endl;
    // for (size_t i = 0; i < std::min(size_t(10), allVertices.size()); i++) {
    //     std::cout << allVertices[i].x << ", " << allVertices[i].y << ", " << allVertices[i].z << std::endl;
    // }

    // std::cout << "First 10 triangles:" << std::endl;
    // for (size_t i = 0; i < std::min(size_t(10), allTriangleIndices.size()); i++) {
    //     std::cout << allTriangleIndices[i].v1 << ", " << allTriangleIndices[i].v2 << ", " << allTriangleIndices[i].v3 << std::endl;
    // }

    // std::vector<Cluster> clusters = createSearchSpaces(allTriangles, weights, tau_S);
    // std::cout << "Finished createSearchSpaces." << std::endl;

    // for (const auto& cluster : clusters) {
    //     if (cluster.triangles.empty()) {
    //         std::cerr << "Found an empty cluster!" << std::endl;
    //     }
    //     for (const auto& triangle : cluster.triangles) {
    //         if (!triangle.isValidity()) {
    //             std::cerr << "Found an invalid triangle!" << std::endl;
    //         }
    //     }
    // }

    // // 2. Randomized Growth Optimization
    // std::cout << "Starting randomizedGrowthOptimization..." << std::endl;
    // std::vector<Component> allComponents;
    // for (const auto& cluster : clusters) {
    //     std::vector<Component> components = randomizedGrowthOptimization(cluster);
    //     allComponents.insert(allComponents.end(), components.begin(), components.end());
    // }
    // std::cout << "Finished randomizedGrowthOptimization." << std::endl;


    // 4. Color Code Search Space Segments (Clusters)
    // float hueStep = 360.0f / clusters.size();  // Divide the hue range (0-360) by the number of clusters
    // float currentHue = 0.0f;
    // float saturation = 1.0f;  // Full saturation for bright colors
    // float value = 1.0f;       // Full value for bright colors

    // for (auto& cluster : clusters) {
    //     Vector3D assignedColor = HSVtoRGB(currentHue, saturation, value);
    //     for (auto& triangle : cluster.triangles) {  // Note: Changed const auto& to auto& since we're modifying the triangle
    //         triangle.setColor(assignedColor);
    //     }
    //     currentHue += hueStep;  // Move to the next hue value
    // }
    // std::cout << "Finished assigning colors to clusters." << std::endl;

    // // 5. Render or save the model with the assigned colors for search space segments
    // std::cout << "Starting saving to OBJ..." << std::endl;
    // saveToOBJ(clusters, "output_filename_without_extension");
    // std::cout << "Finished saving to OBJ." << std::endl;

    std::cout << "Program completed successfully." << std::endl;

    return 0;
}