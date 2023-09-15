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
    std::unordered_map<Triangle, std::array<std::tuple<int, int, int>, 3>> precomputedHashes;
    std::unordered_map<Triangle, std::vector<Triangle>> precomputedNeighbors;

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

    void insert(const Triangle& t, int clusterIndex) {
        auto triangleHashes = getTriangleHashes(t);
        for (const auto& h : triangleHashes) {
            hashTable[h].emplace_back(IndexedTriangle{clusterIndex, t});
        }
        precomputedHashes[t] = triangleHashes;
    }


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

    std::unordered_set<int> getNeighboringClusters(const Triangle& t) {
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
                    clusterIndices.insert(indexedTriangle.clusterIndex);
                }
            }
        }
        return clusterIndices;
    }

    std::vector<Triangle> getPotentialNeighbors(const Triangle& t) {
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
                if (potentialNeighbor != t) {  // Exclude the triangle itself
                    neighbors.insert(potentialNeighbor);
                }
            }
        }
        return std::vector<Triangle>(neighbors.begin(), neighbors.end());
    }


    std::unordered_set<int> getNeighboringClustersForCluster(const Cluster& cluster) {
        std::unordered_set<int> clusterIndices;
        std::unordered_set<Triangle> processedTriangles;
            
        for (const auto& triangle : cluster.triangles) {
            auto neighboringClustersForTriangle = getNeighboringClusters(triangle);
            // std::cout << "For triangle " << triangle.toString() << ", found " << neighboringClustersForTriangle.size() << " neighboring clusters." << std::endl;
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

    void precomputeTriangleNeighbors(const Triangle& triangle) {
        precomputedNeighbors[triangle] = getPotentialNeighbors(triangle);
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

    // Print adjacency list
    // std::cout << "[Debug] Adjacency list:" << std::endl;
    // for (int i = 0; i < adjList.size(); ++i) {
    //     std::cout << i << ": ";
    //     for (int j : adjList[i]) {
    //         std::cout << j << " ";
    //     }
    //     std::cout << std::endl;
    // }

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



struct Component {
    std::vector<Triangle> triangles;
    double weight;

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

};

class Hull {
private:
    std::vector<Point_3> points;
    Polyhedron P;

public:

    Hull(const Component& component) {
        for (const Triangle& triangle : component.triangles) {
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
        // Debugging: Print the number of points.
        //ENABLE DEBUG
        // std::cout << "Number of points before computing hull: " << points.size() << std::endl;
        
        if (points.size() < 4) {
            //ENABLE DEBUG
            // std::cout << "Insufficient points for 3D hull." << std::endl;
            return;
        }
        try {
            CGAL::convex_hull_3(points.begin(), points.end(), P);
        } catch (const std::exception& e) {
            std::cerr << "Exception: " << e.what() << std::endl;
        }

        // Debugging: Print the number of vertices in the computed hull.
        //ENABLE DEBUG
        // std::cout << "Number of vertices in the hull after computing: " << P.size_of_vertices() << std::endl;
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

bool Component::overlapsSignificantly(const Component& other) {
        Hull thisHull = Hull(*this);
        Hull otherHull = Hull(other);

        double thisVolume = thisHull.volume();
        double otherVolume = otherHull.volume();

        Hull combinedHull(thisHull, otherHull);  // You will need to implement this constructor
        double combinedVolume = combinedHull.volume();

        // Compute the volume that is exclusive to each hull
        double exclusiveThisVolume = combinedVolume - otherVolume;
        double exclusiveOtherVolume = combinedVolume - thisVolume;

        // Check for significant overlap
        if (exclusiveThisVolume / thisVolume <= 0.1 || exclusiveOtherVolume / otherVolume <= 0.1) {
            return true;
        }
        return false;

    }

double ATaTb(const Cluster& Ta, const Cluster& Tb, 
             const std::unordered_map<Triangle, std::vector<Triangle>>& potentialNeighborsCache, 
             const TriangleAdjacency& triangleAdjacency,
             const std::unordered_set<Triangle>& tbTriangles) {  // <-- Added this
    int count = 0;

    // No longer converting Tb's triangles to a set here since we're passing it as a parameter

    for (const auto& ta_i : Ta.triangles) {
        auto potentialNeighbors = potentialNeighborsCache.at(ta_i);
        for (const auto& potentialNeighbor : potentialNeighbors) {
            if (tbTriangles.count(potentialNeighbor) && ta_i.isAdjacent(potentialNeighbor, triangleAdjacency)) {
                count++;
            }
        }
    }

    double result = static_cast<double>(count) / std::min(Ta.triangles.size(), Tb.triangles.size());
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

// std::vector<Cluster> initialClusteringByShapeSimilarity(const std::vector<Triangle>& triangles, const std::vector<double>& weights, double tau_S) {
//     // Sort triangles in order of increasing area
//     std::vector<Triangle> sortedTriangles = triangles;
//     std::sort(sortedTriangles.begin(), sortedTriangles.end(), [](const Triangle& a, const Triangle& b) {
//         return a.area < b.area;
//     });

//     std::vector<Cluster> clusters;

//     for (size_t i = 0; i < sortedTriangles.size(); ++i) {
//         bool assigned = false;
//         for (Cluster& cluster : clusters) {
//             double similarity = Stitj(sortedTriangles[i], cluster.triangles[0], weights);
//             if (similarity <= tau_S) {
//                 cluster.triangles.push_back(sortedTriangles[i]);
//                 assigned = true;
//                 break;
//             }
//         }
//         if (!assigned) {
//             // Create a new cluster for this triangle
//             clusters.push_back({{sortedTriangles[i]}});
//         }
//     }
//     return clusters;
// }

std::vector<Cluster> initialClusteringByShapeSimilarity(const std::vector<Triangle>& triangles, const std::vector<double>& weights, double tau_S) {
    // Sort triangles in order of increasing area
    std::vector<Triangle> sortedTriangles = triangles;
    std::sort(sortedTriangles.begin(), sortedTriangles.end(), [](const Triangle& a, const Triangle& b) {
        return a.area < b.area;
    });
    std::vector<double> all_stitj;  // Vector to hold all Stitj values

    std::vector<Cluster> clusters;

    for (size_t i = 0; i < sortedTriangles.size(); ++i) {
        bool assigned = false;
        for (Cluster& cluster : clusters) {
            for (const Triangle& existingTriangle : cluster.triangles) {  // Compare with every triangle in the cluster
                double similarity = Stitj(sortedTriangles[i], existingTriangle, weights);
                all_stitj.push_back(similarity); 
                if (similarity <= tau_S) {
                    cluster.triangles.push_back(sortedTriangles[i]);
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
            clusters.push_back({{sortedTriangles[i]}});
        }
    }
    writeVectorToFile(all_stitj, "all_stitj_non_normalize.txt"); 
    return clusters;
}

bool attemptToMergeComponents(Component& component1, Component& component2, TriangleAdjacency& adjacency) {
    // Create a temporary component with the triangles from both components
    Component tempComponent = component1;
    tempComponent.triangles.insert(tempComponent.triangles.end(), component2.triangles.begin(), component2.triangles.end());

    // Check if the merged component is valid
    if (tempComponent.isValid()) {
        component1 = tempComponent;  // Update component1 to be the merged component
        return true;
    }
    return false;
}

std::vector<std::pair<Triangle, double>> chooseNextTriangleByHullForEachComponent(
    const std::vector<Component>& tempComponents, 
    const std::vector<Triangle>& remainingTriangles
) {
    std::vector<std::pair<Triangle, double>> bestTrianglesForComponents;
    std::unordered_set<Triangle> alreadyChosenTriangles;

    // std::cout << "Number of tempComponents: " << tempComponents.size() << std::endl;
    // std::cout << "Number of remainingTriangles: " << remainingTriangles.size() << std::endl;

    int componentIndex = 0;

    for (const Component& component : tempComponents) {
        Triangle bestTriangle;
        double bestScore = std::numeric_limits<double>::max();
        Hull originalHull(component);
        double originalVolume = originalHull.volume();

        // std::cout << "Iterating for Component Index: " << componentIndex << std::endl;
        // std::cout << "original component" << component.toString() << std::endl;

        for (const Triangle& triangle : remainingTriangles) {
            if (alreadyChosenTriangles.find(triangle) != alreadyChosenTriangles.end()) {
                // std::cout << "Skipping triangle as it's already chosen: " << triangle.toString() << std::endl;
                continue;
            }

            Hull tempHull = originalHull;
            tempHull.add(triangle);
            tempHull.computeHull();
            double newVolume = tempHull.volume();
            double volumeChange = std::abs(newVolume - originalVolume);

            if (volumeChange < bestScore) {
                bestScore = volumeChange;
                bestTriangle = triangle;
            }
        }

        // std::cout << "Component Index: " << componentIndex << ", Best triangle: " 
        //           << bestTriangle.toString() << " with score " << bestScore << std::endl;

        bestTrianglesForComponents.push_back({bestTriangle, bestScore});
        alreadyChosenTriangles.insert(bestTriangle);

        ++componentIndex;
    }

    return bestTrianglesForComponents;
}


std::vector<std::vector<std::pair<Triangle, double>>> getCandidateTrianglesForEachComponent(
    const std::vector<Component>& tempComponents, 
    const std::vector<Triangle>& remainingTriangles,
    size_t topN = 3  // The number of top candidates to consider for each component
) {
    std::vector<std::vector<std::pair<Triangle, double>>> candidateTrianglesForComponents;
    std::unordered_set<Triangle> alreadyChosenTriangles;

    for (const Component& component : tempComponents) {
        std::priority_queue<std::pair<double, Triangle>> candidatePriorityQueue;
        Hull originalHull(component);
        double originalVolume = originalHull.volume();

        for (const Triangle& triangle : remainingTriangles) {
            if (alreadyChosenTriangles.find(triangle) != alreadyChosenTriangles.end()) {
                continue;
            }

            Hull tempHull = originalHull;
            tempHull.add(triangle);
            tempHull.computeHull();
            double newVolume = tempHull.volume();
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


bool recursiveBacktracking(std::vector<Component>& tempComponents,
                           std::unordered_set<Triangle>& remainingTriangles,
                           TriangleAdjacency& adjacency,
                           int depth = 0,
                           const int MAX_DEPTH = 3) {
    std::cout << "Entering recursiveBacktracking at depth: " << depth << std::endl;
    std::cout << "Number of remaining triangles: " << remainingTriangles.size() << std::endl;

    if (depth > MAX_DEPTH) {
        std::cout << "Max depth reached. Returning false." << std::endl;
        return false;
    }

    if (remainingTriangles.empty()) {
        std::cout << "No remaining triangles. Returning true." << std::endl;
        return true;
    }

    auto candidateTrianglesForComponents = getCandidateTrianglesForEachComponent(
        tempComponents,
        std::vector<Triangle>(remainingTriangles.begin(), remainingTriangles.end())
    );

    for (size_t i = 0; i < tempComponents.size(); ++i) {
        auto& component = tempComponents[i];
        Component prevState = component;

        // Retrieve candidate triangles for the current component
        auto& candidateTrianglesForComponent = candidateTrianglesForComponents[i];

        std::cout << "Trying candidate triangles for component index: " << i << std::endl;

        // Loop through the candidate triangles for the current component
        for (const auto& scoreTrianglePair : candidateTrianglesForComponent) {
            Triangle candidateTriangle = scoreTrianglePair.first;

            std::cout << "Trying to add triangle: " << candidateTriangle.toString() << std::endl;
            
            // Try adding each candidate triangle
            if (remainingTriangles.erase(candidateTriangle) == 0) {
                continue;
            }
            component.triangles.push_back(candidateTriangle);

            if (component.isValid()) {
                std::cout << "Component is valid. Recurring..." << std::endl;
                
                if (recursiveBacktracking(tempComponents, remainingTriangles, adjacency, depth + 1)) {
                    return true;
                } else {
                    std::cout << "Backtracking due to unsuccessful recursion..." << std::endl;
                }
            } else {
                std::cout << "Component is not valid. Rolling back..." << std::endl;
            }

            // If you reach here, it means you need to backtrack
            component = prevState;
            remainingTriangles.insert(candidateTriangle);
        }
    }

    std::cout << "All paths explored. Returning false." << std::endl;
    return false;
}





std::vector<Triangle> findLargestDisjointSimilarTriangles(const std::vector<Triangle>& triangles, std::vector<double>& weights, double tau_S) {
    std::vector<Triangle> seeds;

    // Sort triangles by area
    std::vector<Triangle> sortedTriangles = triangles;
    std::sort(sortedTriangles.begin(), sortedTriangles.end(), [](const Triangle& a, const Triangle& b) {
        return a.area > b.area;
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

    // Find and print out the largest seed
    auto largestSeedIter = std::max_element(seeds.begin(), seeds.end(), [](const Triangle& a, const Triangle& b) {
        return a.area < b.area;
    });

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
            Component newComponent;
            newComponent.triangles.push_back(seed);
            tempComponents.push_back(newComponent);
        }
    }
    return tempComponents;
}

bool growComponents(std::vector<Component>& tempComponents,
                    std::unordered_set<Triangle>& remainingTriangles,
                    TriangleAdjacency& adjacency) {
    bool anyComponentGrew = false;

    // std::cout << "[Debug] Starting growComponents. remainingTriangles size: " << remainingTriangles.size() << std::endl;

    if(remainingTriangles.empty()) {
        // std::cout << "[Debug] remainingTriangles is already empty at the start. Exiting.\n";
        return anyComponentGrew; // You could return false or an appropriate value
    }

    auto bestTrianglesForComponents = chooseNextTriangleByHullForEachComponent(
        tempComponents, 
        std::vector<Triangle>(remainingTriangles.begin(), remainingTriangles.end())
    );

    for (size_t i = 0; i < tempComponents.size(); ++i) {
        Component& component = tempComponents[i];
        Triangle nextTriangle = std::get<0>(bestTrianglesForComponents[i]);

        if (remainingTriangles.erase(nextTriangle) == 0) {
            // std::cout << "[Debug] Triangle not found in remainingTriangles, skipping" << std::endl;
            continue;
        }

        // std::cout << "[Debug] Removed triangle, new remainingTriangles size: " << remainingTriangles.size() << std::endl;

        // Keep a snapshot of the component before attempting to grow it
        Component prevState = component;

        if (component.isValid()) {
            // std::cout << "[Debug] Valid growth BEFORE PUSHING NEW TRIANGLE.\n";
        } 

        // Attempt to grow the component
        component.triangles.push_back(nextTriangle);

        if (component.isValid()) {
            // std::cout << "[Debug] Valid growth.\n";
            anyComponentGrew = true;
        } else {
            // std::cout << "[Debug] Invalid growth. Reverting...\n";
            // If the growth was invalid, revert to the previous state
            component = prevState;
            remainingTriangles.insert(nextTriangle);
            // std::cout << "[Debug] Reverted. New remainingTriangles size: " << remainingTriangles.size() << std::endl;
        }
    }

    // std::cout << "[Debug] Ending growComponents. remainingTriangles size: " << remainingTriangles.size() << std::endl;

    return anyComponentGrew;
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

std::vector<Component> randomizedGrowthOptimization(const Cluster& cluster,
                                                    TriangleAdjacency& adjacency,
                                                    std::vector<double>& weights,
                                                    double tau_S) {
    static int iteration = 1;  // Keeping it static as in your original code
    std::vector<Component> finalComponents;  // To hold the final components
    std::unordered_set<Triangle> remainingTriangles(cluster.triangles.begin(), cluster.triangles.end());

    while (!remainingTriangles.empty()) {
        // std::cout << "remaining triangles size: " << remainingTriangles.size() << std::endl;
        auto start = std::chrono::high_resolution_clock::now();  // Start timing

        std::vector<Component> currentIterationComponents;  // Components for this iteration

        // Seed with largest disjoint similar triangles
        std::vector<Triangle> seeds = findSeeds(remainingTriangles, weights, tau_S);
        std::vector<Component> tempComponents = initializeTempComponents(seeds, remainingTriangles);

        // Sort tempComponents based on the updated weights
        std::sort(tempComponents.begin(), tempComponents.end(), [](const Component& a, const Component& b) {
            double weight_a = 1.0 / (a.triangles.size() * a.triangles.size());
            double weight_b = 1.0 / (b.triangles.size() * b.triangles.size());
            return weight_a < weight_b;
        });

        bool canGrow = true;

        while (canGrow) {
            canGrow = growComponents(tempComponents, remainingTriangles, adjacency);

            // After growing components and before the next iteration
            // Perform the merging of significantly overlapping components
            for (int i = 0; i < tempComponents.size(); ++i) {
                for (int j = i + 1; j < tempComponents.size(); ++j) {
                    if (tempComponents[i].overlapsSignificantly(tempComponents[j])) {
                        // Attempt to merge the components (implement this based on your own criteria)
                        // Also remove one of them from the tempComponents list if merge is successful
                        if (attemptToMergeComponents(tempComponents[i], tempComponents[j], adjacency)) {
                            tempComponents.erase(tempComponents.begin() + j);
                            --j;
                        }
                    }
                }
            }

            // Evaluate components based on the paper's criteria
            for (auto& component : tempComponents) {
                component.updateWeight();  // Assuming you have a method that updates the weight
            }

            // Sort based on weight
            std::sort(tempComponents.begin(), tempComponents.end(), [](const Component& a, const Component& b) {
                return a.weight < b.weight;
            });

            // Check contiguity of each component
            for (const auto& component : tempComponents) {
                if (!areTrianglesContiguous(component.triangles)) {
                    std::cout << "Component became non-contiguous during growth." << std::endl;
                }
            }

            if (!canGrow) {
                // std::cout << "remaining triangles size before backtracing: " << remainingTriangles.size() << std::endl;
                if (!recursiveBacktracking(tempComponents, remainingTriangles, adjacency)) {
                    // No more alternative paths to explore; break out of loop
                    // std::cout << "Backtracking failed. Breaking out." << std::endl;
                    break;
                }
            }
        }

        auto end = std::chrono::high_resolution_clock::now();  // End timing
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "Time taken by loop iteration: " << duration.count() << " microseconds" << std::endl;

        // Validate and add valid components to this iteration's components
        for (const auto& comp : tempComponents) {
            if (comp.isValid()) {
                currentIterationComponents.push_back(comp);
            }
        }

        // Add this iteration's components to the final list
        finalComponents.insert(finalComponents.end(),
                               currentIterationComponents.begin(),
                               currentIterationComponents.end());

        std::cout << "Iteration " << iteration << " completed." << std::endl;
        iteration++;

        // If there are leftover triangles, the loop will automatically re-seed in the next iteration
    }

    // Final validation
    if (!validateFinalState(finalComponents, cluster)) {
        std::cout << "Validation failed." << std::endl;
    }

    return finalComponents;
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

std::unordered_map<Triangle, size_t> triangleToClusterMap;  // Global or member variable

// Fast lookup function
size_t findClusterIndexForTriangle(const Triangle& targetTriangle) {
    auto it = triangleToClusterMap.find(targetTriangle);
    if (it != triangleToClusterMap.end()) {
        return it->second;
    }
    return static_cast<size_t>(-1); // Not found
}

// Update the map whenever you move a triangle from one cluster to another
void moveTriangleToCluster(const Triangle& triangle, size_t newClusterIndex) {
    triangleToClusterMap[triangle] = newClusterIndex;
}

void checkForDuplicates(const std::vector<Cluster>& clusters) {
    for (const auto& cluster : clusters) {
        std::unordered_set<Triangle> uniqueTriangles(cluster.triangles.begin(), cluster.triangles.end());
        if (uniqueTriangles.size() != cluster.triangles.size()) {
            std::cerr << "Duplicate triangles found in a cluster" << std::endl;
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

bool checkConsistency(
    const std::vector<Cluster>& clusters, 
    size_t expectedTotalTriangles, 
    const std::string& operation, 
    size_t sourceClusterIndex, 
    size_t targetClusterIndex,
    size_t trianglesAdded,
    size_t trianglesRemoved
) {
    size_t actualTotalTriangles = 0;
    for (size_t i = 0; i < clusters.size(); ++i) {
        const auto& cluster = clusters[i];
        actualTotalTriangles += cluster.triangles.size();
    }

    if (actualTotalTriangles != expectedTotalTriangles) {
        std::cerr << "[Operation: " << operation << "] Inconsistency: expected " << expectedTotalTriangles 
                  << ", got " << actualTotalTriangles
                  << ". Source Cluster: " << sourceClusterIndex
                  << ". Target Cluster: " << targetClusterIndex
                  << ". Triangles Added: " << trianglesAdded
                  << ". Triangles Removed: " << trianglesRemoved
                  << std::endl;
        return false;  // Inconsistency found
    }

    return true;  // No inconsistency
}

// Define a threshold for the variance below which you consider the system to have converged
double convergenceMetricThreshold = 0.1;  // Replace with your desired threshold value

void refineClusters(std::vector<Cluster>& clusters, double tau_N) {
    int iteration = 0;
    SpatialHash spatialHash(0.1);
    std::deque<int> lastNMergedAndSplitSums;
    int N = 3;  // Number of last iterations to consider
    int changeRateThreshold = 20;  // The threshold for the difference in the sum of merged and split triangles
    auto start = std::chrono::high_resolution_clock::now();
    spatialHash.precomputeTriangleHashes(clusters);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Time taken for precomputeTriangleHashes: " << elapsed.count() << " seconds" << std::endl;

    size_t initialTriangleCount = 0;
    for (const auto& cluster : clusters) {
        initialTriangleCount += cluster.triangles.size();
    }

    while (true) {
        // 1. Initialize the TriangleAdjacency object
        TriangleAdjacency triangleAdjacency;
        for (const auto& cluster : clusters) {
            for (const auto& triangle : cluster.triangles) {
                triangleAdjacency.addTriangle(triangle);
            }
        }

        std::cout << "Starting iteration " << iteration << " with " << clusters.size() << " clusters." << std::endl;

        int mergesInThisIteration = 0;
        int splitsInThisIteration = 0;
        size_t pairsProcessed = 0;  // Initialize counter for this iteration

        // Declare a variable to accumulate the total time for the specific loop
        std::chrono::milliseconds totalLoopTime(0);

        std::map<std::pair<size_t, size_t>, std::vector<Triangle>> mergeMap;
        std::map<std::pair<size_t, size_t>, std::vector<Triangle>> splitMap;

        std::vector<std::unordered_set<Triangle>> clusterTrianglesSet(clusters.size());
        for (size_t i = 0; i < clusters.size(); ++i) {
            clusterTrianglesSet[i] = std::unordered_set<Triangle>(clusters[i].triangles.begin(), clusters[i].triangles.end());
        }


        for (size_t i = 0; i < clusters.size(); ++i) {

            if (clusters[i].triangles.empty()) continue;

            std::unordered_set<int> neighboringClusterIndices = spatialHash.getNeighboringClustersForCluster(clusters[i]);
            std::unordered_map<Triangle, std::vector<Triangle>> potentialNeighborsCache;

            for (const auto& triangleA : clusters[i].triangles) {
                potentialNeighborsCache[triangleA] = spatialHash.getPotentialNeighbors(triangleA);
            }

            std::vector<int> neighboringClusterIndicesVec(neighboringClusterIndices.begin(), neighboringClusterIndices.end());

            auto loopStartTime = std::chrono::high_resolution_clock::now();  // Start time measurement
            std::unordered_set<Triangle> alreadyAddedToMaps;  // New line


            // #pragma omp parallel for  // Enable OpenMP parallelization
            for (size_t idx = 0; idx < neighboringClusterIndicesVec.size(); ++idx) {
                int j = neighboringClusterIndicesVec[idx];
                if (i == j || clusters[j].triangles.empty()) continue;

                std::unordered_set<Triangle> localMergeList;
                std::unordered_set<Triangle> tbTriangles(clusters[j].triangles.begin(), clusters[j].triangles.end());

                double similarity = ATaTb(clusters[i], clusters[j], potentialNeighborsCache, triangleAdjacency, tbTriangles);
                std::unordered_set<Triangle> toBeMerged;
                std::unordered_set<Triangle> toBeSplit;

                if (similarity >= tau_N) {
                    for (const auto& triangleA : clusters[i].triangles) {

                        bool foundAdjacent = false;
                        for (const auto& triangleB : potentialNeighborsCache[triangleA]) {
                            if (triangleA.isAdjacent(triangleB, triangleAdjacency)) {
                                // Both triangles meet the merge condition, so add them to the list
                                toBeMerged.insert(triangleA);
                                toBeMerged.insert(triangleB);  // Insert triangleB as well
                                foundAdjacent = true;
                                break;
                            }
                        }
                        if (!foundAdjacent) {
                            toBeSplit.insert(triangleA);
                        }
                    }

                    for (const auto& triangleB : clusters[j].triangles) {
                            if (toBeMerged.count(triangleB) == 0) {
                                // Only add to be split if it's not in the merge list
                                toBeSplit.insert(triangleB);
                            }
                    }
                    // #pragma omp critical  // Critical section
                    {
                        if (!toBeMerged.empty()) {
                            for (const auto& triangle : toBeMerged) {
                                if (alreadyAddedToMaps.count(triangle) == 0) {
                                    // Add to mergeMap only if not already added
                                    mergeMap[{i, j}].push_back(triangle);
                                    alreadyAddedToMaps.insert(triangle);
                                }
                            }
                        }
                        
                        if (!toBeSplit.empty()) {
                            for (const auto& triangle : toBeSplit) {
                                if (alreadyAddedToMaps.count(triangle) == 0) {
                                    // Add to splitMap only if not already added
                                    splitMap[{i, j}].push_back(triangle);
                                    alreadyAddedToMaps.insert(triangle);
                                }
                            }
                        }
                    }
                }
            }
            auto loopEndTime = std::chrono::high_resolution_clock::now();  // End time measurement
            totalLoopTime += std::chrono::duration_cast<std::chrono::milliseconds>(loopEndTime - loopStartTime);
        }

        std :: cout << "Time taken for the merge loop: " << totalLoopTime.count() << " milliseconds" << std::endl;

        size_t initialCount = 0;
        for (const auto& cluster : clusters) {
            initialCount += cluster.triangles.size();
        }
        std::cout << "Initial Triangle Count: " << initialCount << std::endl;

        // Initialize a variable to keep track of the total number of triangles after all merges
        int postMergeTotalTriangles = 0;
        std::unordered_set<Triangle> mergedTriangles;  // To keep track of merged triangles
        std::unordered_set<size_t> mergedClusters;
        std::unordered_map<Triangle, size_t> originalClusterOfTriangle;  // Triangle -> Original cluster index

        initializeOriginalClusterOfTriangle(originalClusterOfTriangle, clusters);

        // Create a temporary copy of the clusters
        std::vector<Cluster> tempClusters = clusters;
        clusters.clear();

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
                        targetTriangles.push_back(triangle);
                        targetClusterTriangles.insert(triangle);
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

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distr(0, clusters.size() - 1);
        std::unordered_set<size_t> emptiedClusters;  // Keep track of emptied clusters
        std::unordered_map<Triangle, size_t> triangleMovedToCluster;

        for (auto& entry : splitMap) {
            auto& trianglesToSplit = entry.second;
            size_t originalClusterIndex = originalClusterOfTriangle[trianglesToSplit[0]];

            if (emptiedClusters.count(originalClusterIndex) > 0) {
                continue;
            }

            auto& originalClusterTriangles = clusters[originalClusterIndex].triangles;

            size_t notMovedCount = 0;
            std::vector<Triangle> successfullyMoved;

            for (const auto& triangle : originalClusterTriangles) {
                if (movedTriangles.count(triangle) == 0) {
                    successfullyMoved.push_back(triangle);
                    movedTriangles.insert(triangle);  // Mark as moved, but don't actually move yet
                } else {
                    notMovedCount++;
                }
            }

            emptiedClusters.insert(originalClusterIndex);  // Mark this cluster as to be emptied, but don't actually clear it yet
        }


        // Remove empty clusters after merges and splits
        auto itSplit = std::remove_if(clusters.begin(), clusters.end(),
            [](const Cluster& cluster) { return cluster.triangles.empty(); });
            
        std::cout << "Number of empty clusters: " << std::distance(itSplit, clusters.end()) << std::endl;

        clusters.erase(itSplit, clusters.end());


        // Your function to check for duplicates (Assuming you have this implemented)
        checkForDuplicates(clusters);



        size_t finalTriangleCount = 0;
        for (const auto& cluster : clusters) {
            finalTriangleCount += cluster.triangles.size();
        }

        // Verify
        if (initialTriangleCount != finalTriangleCount || trianglesAdded != trianglesRemoved) {
            std::cerr << "ERROR: Triangle count mismatch. Initial: " << initialTriangleCount << ", Final: " << finalTriangleCount << std::endl;
            std::cerr << "Triangles Added: " << trianglesAdded << ", Triangles Removed: " << trianglesRemoved << std::endl;
        }

        validateTotalTriangles(clusters, initialTriangleCount);

        // std::cout << "Total time taken: " << totalLoopSplitTime.count() << " microseconds" << std::endl;

        std::cout << "After " << iteration << " iterations, merged " << mergeMap.size() << " triangles and split " << splitMap.size() << " triangles." << std::endl;

        // After you calculate the merged and split triangles for this iteration
        int mergedAndSplitThisIteration = mergeMap.size() + splitMap.size();
        std::cout << "Merged and split this iteration: " << mergedAndSplitThisIteration << std::endl;

        // Update the deque
        lastNMergedAndSplitSums.push_back(mergedAndSplitThisIteration);
        if (lastNMergedAndSplitSums.size() > N) {
            lastNMergedAndSplitSums.pop_front();
        }

        // Debug print the deque
        std::cout << "last N merged and Split Sums: ";
        for (const auto& val : lastNMergedAndSplitSums) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "last N merged and Split Sums size: " << lastNMergedAndSplitSums.size() << std::endl;

        if (lastNMergedAndSplitSums.size() >= N) { // make sure there are at least N elements
            int sumOfAbsDifferences = 0;
            
            // Loop starts from the Nth element from the end and goes to the last element
            auto itEnd = std::prev(lastNMergedAndSplitSums.end());
            for (auto it = std::prev(lastNMergedAndSplitSums.end(), N); it != itEnd; ++it) {
                // Debug print the current and next elements
                std::cout << "Current element: " << *it << ", Next element: " << *std::next(it) << std::endl;
                
                // Calculate the absolute difference between the current and next element
                int absDifference = std::abs(*it - *std::next(it));
                
                // Debug print the calculated absolute difference
                std::cout << "Calculated absolute difference: " << absDifference << std::endl;
                
                // Add the absolute difference to the sum
                sumOfAbsDifferences += absDifference;
            }
            std :: cout << "sum of abs differences: " << sumOfAbsDifferences << std::endl;
            // Now sumOfAbsDifferences contains the sum of the absolute differences of the last N elements
            if (sumOfAbsDifferences <= changeRateThreshold) {
                std::cout << "Convergence achieved based on change rate. Stopping..." << std::endl;
                break;
            }
            else {
                std::cout << "convergence not achieved based on change rate. Continuing..." << std::endl;
            }
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


std::vector<Cluster> createSearchSpaces(const std::vector<Triangle>& triangles, const std::vector<double>& weights, double tau_S) {    // 1. Initial clustering based on shape similarity
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
    refineClusters(clusters, tau_N);

    std::cout << "Number of clusters after refining: " << clusters.size() << std::endl;

    // std::set<Triangle> usedTriangles; // This will store triangles that are part of a component
    std::vector<Component> allComponents; // This will store all grown components
    
    std::vector<double> weightsPostClustering = {.25, .25, .25, 0};

    // // 4. Randomized Growth Optimization
    std::cout << "4. Randomized Growth Optimization" << std::endl;
    TriangleAdjacency adjacencyInstance; // Assume you've created and possibly initialized this instance
    for (const Cluster& cluster : clusters) {
        for (const Triangle& triangle : cluster.triangles) {
            adjacencyInstance.addTriangle(triangle);
        }
    }

    for (Cluster& cluster : clusters) {
        std::vector<Component> components = randomizedGrowthOptimization(cluster, adjacencyInstance, weightsPostClustering, tau_S);
        allComponents.insert(allComponents.end(), components.begin(), components.end());
    }

    std::cout << "Total components after randomized growth: " << allComponents.size() << std::endl;


    std::cout << "Before assignContrastingColors" << std::endl;
    assignContrastingColors(clusters);
    std::cout << "After assignContrastingColors" << std::endl;

    std::cout << "Before assignContrastingColorsToComponents" << std::endl;
    assignContrastingColorsToComponents(allComponents);
    std::cout << "After assignContrastingColorsToComponents" << std::endl;

    std::cout << "Before saveToOBJ" << std::endl;
    saveToOBJ(clusters, "clusters");
    std::cout << "After saveToOBJ" << std::endl;

    std::cout << "Before saveComponentsToOBJ" << std::endl;
    saveComponentsToOBJ(allComponents, "components");
    std::cout << "After saveComponentsToOBJ" << std::endl;

    std::cout << "Saved refined clusters to OBJ for visualization." << std::endl;

}



void testSpatialHash(const std::vector<Triangle>& triangles) {

    SpatialHash spatialHash(0.1);  // Assuming 1 is the cell size you want

    // Ensure there are at least 3 triangles for the test
    if (triangles.size() < 3) {
        std::cout << "Not enough triangles for the test." << std::endl;
        return;
    }

    // Use the first 3 triangles for the test
    std::cout << "Before triangle 1 assignment" << std::endl;
    Triangle t1 = triangles[0];
    std::cout << "After triangle 1 assignment" << std::endl;
    Triangle t2 = triangles[1];
    Triangle t3 = triangles[2];

    spatialHash.insert(t1, 0);
    spatialHash.insert(t2, 1);
    spatialHash.insert(t3, 2);

    auto neighborsT1 = spatialHash.getPotentialNeighbors(t1);
    auto neighborsT2 = spatialHash.getPotentialNeighbors(t2);
    auto neighborsT3 = spatialHash.getPotentialNeighbors(t3);

    std::cout << "Neighbors for T1: " << neighborsT1.size() << std::endl;
    std::cout << "Neighbors for T2: " << neighborsT2.size() << std::endl;
    std::cout << "Neighbors for T3: " << neighborsT3.size() << std::endl;
}

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

int main() {
    std::cout << "Starting program..." << std::endl;
    // std::vector<Triangle> allTriangles = {};  // Fill with triangles from your data source
    std::vector<double> weights = {.3, .3, .3, 0};
    std::vector<Vector3D> allVertices;    
    double tau_S =.01;  // Threshold for shape similarity
    std::string inputfile = "/home/kingpin/Documents/blender/burger.obj";
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
    std::vector<double> maxEdges;
    std::vector<double> minEdges;
    for(const Triangle& t : allTriangles) {
        areas.push_back(t.area);
        maxEdges.push_back(t.e_max);
        minEdges.push_back(t.e_min);
    }

    std::sort(areas.begin(), areas.end());

    double minValue = *std::min_element(areas.begin(), areas.end());
    double maxValue = *std::max_element(areas.begin(), areas.end());

    double minValueEdgeMin = *std::min_element(minEdges.begin(), minEdges.end());
    double maxValueEdgeMin = *std::max_element(minEdges.begin(), minEdges.end());

    double minValueEdgeMax = *std::min_element(maxEdges.begin(), maxEdges.end());
    double maxValueEdgeMax = *std::max_element(maxEdges.begin(), maxEdges.end());

    // // Compute statistics
    auto computeStats = [](const std::vector<double>& values) {
        double min_val = *std::min_element(values.begin(), values.end());
        double max_val = *std::max_element(values.begin(), values.end());
        double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
        double variance = std::accumulate(values.begin(), values.end(), 0.0, [mean](double acc, double val) {
            return acc + std::pow(val - mean, 2);
        }) / values.size();
        double std_dev = std::sqrt(variance);
        
        return std::make_tuple(min_val, max_val, mean, std_dev);
    };
    
    auto [min_area, max_area, mean_area, std_dev_area] = computeStats(areas);
    auto [min_max_edge, max_max_edge, mean_max_edge, std_dev_max_edge] = computeStats(maxEdges);
    auto [min_min_edge, max_min_edge, mean_min_edge, std_dev_min_edge] = computeStats(minEdges);
    
    // // Print statistics
    // std::cout << "Area: Min=" << min_area << ", Max=" << max_area << ", Mean=" << mean_area << ", StdDev=" << std_dev_area << std::endl;
    // std::cout << "Max Edge: Min=" << min_max_edge << ", Max=" << max_max_edge << ", Mean=" << mean_max_edge << ", StdDev=" << std_dev_max_edge << std::endl;
    // std::cout << "Min Edge: Min=" << min_min_edge << ", Max=" << max_min_edge << ", Mean=" << mean_min_edge << ", StdDev=" << std_dev_min_edge << std::endl;

    for (auto& triangle : allTriangles) {
        triangle.area = minMaxScale(triangle.area, min_area, max_area);
        triangle.e_max = minMaxScale(triangle.e_max, min_max_edge, max_max_edge);
        triangle.e_min = minMaxScale(triangle.e_min, min_min_edge, max_min_edge);
    }

    areas.clear();
    maxEdges.clear();
    minEdges.clear();
    for (const Triangle& t : allTriangles) {
        areas.push_back(t.area);
        maxEdges.push_back(t.e_max);
        minEdges.push_back(t.e_min);
    }

    auto [min_area_norm, max_area_norm, mean_area_norm, std_dev_area_norm] = computeStats(areas);
    auto [min_max_edge_norm, max_max_edge_norm, mean_max_edge_norm, std_dev_max_edge_norm] = computeStats(maxEdges);
    auto [min_min_edge_norm, max_min_edge_norm, mean_min_edge_norm, std_dev_min_edge_norm] = computeStats(minEdges);

    // std::cout << "After normalization:" << std::endl;
    std::cout << "Area: Min=" << min_area_norm << ", Max=" << max_area_norm << ", Mean=" << mean_area_norm << ", StdDev=" << std_dev_area_norm << std::endl;
    std::cout << "Max Edge: Min=" << min_max_edge_norm << ", Max=" << max_max_edge_norm << ", Mean=" << mean_max_edge_norm << ", StdDev=" << std_dev_max_edge_norm << std::endl;
    std::cout << "Min Edge: Min=" << min_min_edge_norm << ", Max=" << max_min_edge_norm << ", Mean=" << mean_min_edge_norm << ", StdDev=" << std_dev_min_edge_norm << std::endl;
    

    std::cout << "First 10 triangles:" << std::endl;
    for (int i = 0; i < 10 && i < allTriangles.size(); i++) {
        std::cout << "Triangle " << i << ": v1(" << allTriangles[i].vec1.x << "," << allTriangles[i].vec1.y << "," << allTriangles[i].vec1.z 
                << ") v2(" << allTriangles[i].vec2.x << "," << allTriangles[i].vec2.y << "," << allTriangles[i].vec2.z 
                << ") v3(" << allTriangles[i].vec3.x << "," << allTriangles[i].vec3.y << "," << allTriangles[i].vec3.z << ")" << std::endl;
    }

    testSpatialHash(allTriangles);

    std::cout << "Number of triangles to be clustered: " << allTriangles.size() << std::endl;

    std::vector<Cluster> clusters = createSearchSpaces(allTriangles, weights, tau_S); 
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