#include <vector>
#include <algorithm>
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
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
// #include "/tinyobj_loader_opt.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Polyhedron_3<K> Polyhedron;
class Triangle;
class TriangleAdjacency;
class Component; // forward declare the Component class
// bool areTrianglesContiguous(const std::vector<Triangle>& triangles);
// bool areComponentsSimilar(const Component& a, const Component& b, double similarityThreshold);
// void saveTrianglesToOBJ(const std::vector<Triangle>& triangles, const std::string& filename);

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
};

// Triangle methods
Triangle::Triangle(const Vector3D& a, const Vector3D& b, const Vector3D& c)
    : vec1(a), vec2(b), vec3(c) {
    Vector3D edge1 = vec2 - vec1;
    Vector3D edge2 = vec3 - vec1;
    normal = edge1.cross(edge2);
    double normalLength = normal.length();
    normal = normal / normalLength;
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





// class TriangleAdjacency {
// private:
//     typedef TriangleEdge Edge;
//     std::unordered_map<Edge, std::vector<Triangle>, EdgeHash> adjacencyMap;

// public:
//     void addTriangle(const Triangle& triangle);
//     bool isAdjacent(const Triangle& t1, const Triangle& t2) const;
// };

// class Triangle {
// public:
//     Vector3D vec1, vec2, vec3;  // Geometric vertices
//     Vector3D normal;
//     double area;
//     double e_min;
//     double e_max;
//     Vector3D color;  // RGB color for the triangle

//     // Default constructor
//     Triangle() : area(0), e_min(0), e_max(0) {}

//     // Constructor using Vector3D
//     Triangle(const Vector3D& a, const Vector3D& b, const Vector3D& c)
//         : vec1(a), vec2(b), vec3(c) {
//         // Compute edges
//         Vector3D edge1 = vec2 - vec1;
//         Vector3D edge2 = vec3 - vec1;

//         // Compute normal
//         normal = edge1.cross(edge2);
//         double normalLength = normal.length();
//         normal = normal / normalLength;  // Normalize the normal vector

//         // Compute area
//         area = 0.5 * normalLength;

//         // Compute edge lengths
//         double len1 = edge1.length();
//         double len2 = edge2.length();
//         double len3 = (vec3 - vec2).length();
//         e_min = std::min({len1, len2, len3});
//         e_max = std::max({len1, len2, len3});
//     }

//     // Check if triangle is valid based on vertices
//     bool isValid() const {
//         return !vec1.isZero() || !vec2.isZero() || !vec3.isZero();
//     }

//     bool operator==(const Triangle& other) const {
//         return (vec1 == other.vec1 && vec2 == other.vec2 && vec3 == other.vec3) ||
//                (vec1 == other.vec1 && vec2 == other.vec3 && vec3 == other.vec2) ||
//                (vec1 == other.vec2 && vec2 == other.vec1 && vec3 == other.vec3) ||
//                (vec1 == other.vec2 && vec2 == other.vec3 && vec3 == other.vec1) ||
//                (vec1 == other.vec3 && vec2 == other.vec1 && vec3 == other.vec2) ||
//                (vec1 == other.vec3 && vec2 == other.vec2 && vec3 == other.vec1);
//     }

//     bool operator!=(const Triangle& other) const {
//         return !(*this == other);
//     }

//     bool operator<(const Triangle& other) const {
//         return this->area < other.area;
//     }

//     void setColor(const Vector3D& assignedColor) {
//         color = assignedColor;
//     }

//     std::vector<TriangleEdge> getEdges() const {
//         return {TriangleEdge(vec1, vec2), TriangleEdge(vec2, vec3), TriangleEdge(vec3, vec1)};
//     }

//     bool isAdjacent(const Triangle& other, const TriangleAdjacency* adjacency) const {
//         if (adjacency) {
//             return adjacency->isAdjacent(*this, other);
//         }
//         return false; // or throw an error if adjacency is nullptr
//     }


// };


// class TriangleAdjacency {
// private:
//     typedef TriangleEdge Edge;
//     std::unordered_map<Edge, std::vector<Triangle>, EdgeHash> adjacencyMap;

// public:
//     void addTriangle(const Triangle& triangle) {
//         // Extract edges from triangle
//         std::vector<TriangleEdge> edges = {
//             TriangleEdge(triangle.vec1, triangle.vec2),
//             TriangleEdge(triangle.vec2, triangle.vec3),
//             TriangleEdge(triangle.vec3, triangle.vec1)
//         };

//         for (const TriangleEdge& edge : edges) {
//             adjacencyMap[edge].push_back(triangle);
//         }
//     }

//     bool isAdjacent(const Triangle& t1, const Triangle& t2) const {
//         std::vector<TriangleEdge> edges = {
//             TriangleEdge(t1.vec1, t1.vec2),
//             TriangleEdge(t1.vec2, t1.vec3),
//             TriangleEdge(t1.vec3, t1.vec1)
//         };

//         for (const TriangleEdge& edge : edges) {
//             if (adjacencyMap.count(edge)) {
//                 const auto& triangles = adjacencyMap.at(edge);
//                 if (std::find(triangles.begin(), triangles.end(), t2) != triangles.end()) {
//                     return true;
//                 }
//             }
//         }

//         return false;
//     }
// };

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

struct MergeSplitDecision {
    size_t clusterA;
    size_t clusterB;
    Cluster mergedCluster;
    Cluster splitClusterA;
    Cluster splitClusterB;
};

std::vector<MergeSplitDecision> decisions;

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
        std::vector<Triangle> neighbors;
            
        // Check if triangle's hashes are precomputed.
        if (precomputedHashes.find(t) == precomputedHashes.end()) {
            std::cout << "Error: Triangle hashes not precomputed for: " << t.toString() << std::endl;
            return neighbors;
        }

        // Retrieve precomputed hash values for the triangle.
        auto hashes = precomputedHashes[t];

        for (const auto& h : hashes) {

            for (const auto& indexedTriangle : hashTable[h]) {
                const Triangle& potentialNeighbor = indexedTriangle.triangle;


                if (potentialNeighbor != t) {  // Exclude the triangle itself
                    neighbors.push_back(potentialNeighbor);
                }
            }
        }
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());

        return neighbors;
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




// std::vector<Triangle> getAdjacentTriangles(const Triangle& target, const std::vector<Triangle>& triangles) {
//     std::vector<Triangle> adjacentTriangles;

//     for (const Triangle& triangle : triangles) {
//         if (areAdjacent(target, triangle)) {  // This function checks if two triangles share an edge or vertex
//             adjacentTriangles.push_back(triangle);
//         }
//     }

//     return adjacentTriangles;
// }

// bool areTrianglesContiguous(const std::vector<Triangle>& triangles) {
//     if (triangles.empty()) return false;

//     std::unordered_set<Triangle> visited;
//     std::stack<Triangle> stack;

//     stack.push(triangles[0]);

//     while (!stack.empty()) {
//         Triangle current = stack.top();
//         stack.pop();

//         if (visited.find(current) == visited.end()) {
//             visited.insert(current);

//             for (Triangle adjacent : getAdjacentTriangles(current, triangles)) {  
//                 if (visited.find(adjacent) == visited.end()) {
//                     stack.push(adjacent);
//                 }
//             }
//         }
//     }

//     return visited.size() == triangles.size();
// }

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

// bool areAdjacent(const Triangle& t1, const Triangle& t2) {
//     if (t1 == t2) {
//         return false;
//     }

//     // Put vertices of t2 in a set for faster lookup
//     std::unordered_set<Vector3D, Vector3DHash> t2Vertices = {t2.vec1, t2.vec2, t2.vec3};

//     // Count the vertices of t1 that are in t2
//     int sharedVerticesCount = 0;
//     sharedVerticesCount += t2Vertices.count(t1.vec1);
//     sharedVerticesCount += t2Vertices.count(t1.vec2);
//     sharedVerticesCount += t2Vertices.count(t1.vec3);

//     // If they share exactly two vertices, they are adjacent
//     return sharedVerticesCount == 2;
// }




// struct Component {
//     std::vector<Triangle> triangles;

//     bool isValid() const {
//         // Check if triangles are contiguous
//         // This requires a function that checks if a group of triangles form a connected set
//         if (!areTrianglesContiguous(triangles)) {  // To be implemented
//             return false;
//         }
//         // Check other properties (iii and iv) if needed
//         return true;
//     }
    
//     double weight() const {
//         return 1.0 / (triangles.size() * triangles.size());
//     }
// };

// class Hull {
// private:
//     std::vector<Point_3> points;
//     Polyhedron P;

// public:

//     Hull(const Component& component) {
//         for (const Triangle& triangle : component.triangles) {
//             add(triangle);
//         }
//     }


//     // Add a triangle to the collection of points
//     void add(const Triangle& triangle) {
//         points.push_back(Point_3(triangle.vec1.x, triangle.vec1.y, triangle.vec1.z));
//         points.push_back(Point_3(triangle.vec2.x, triangle.vec2.y, triangle.vec2.z));
//         points.push_back(Point_3(triangle.vec3.x, triangle.vec3.y, triangle.vec3.z));
//     }

//     // Compute the convex hull and store it in Polyhedron P
//     void computeHull() {
//         CGAL::convex_hull_3(points.begin(), points.end(), P);
//     }

//     double volume() const {
//         double total_volume = 0.0;

//         for (auto facet = P.facets_begin(); facet != P.facets_end(); ++facet) {
//             auto h = facet->halfedge();
            
//             Point_3 A = h->vertex()->point();
//             Point_3 B = h->next()->vertex()->point();
//             Point_3 C = h->opposite()->vertex()->point();

//             // Compute Q_F, we can take A as a point on face F
//             Point_3 Q_F = A;

//             // Compute the normal of the face using the correct vector type
//             K::Vector_3 N_F = CGAL::normal(A, B, C);

//             // Compute area of triangle
//             double areaF = std::sqrt(N_F.squared_length()) * 0.5;
//             N_F = N_F / std::sqrt(N_F.squared_length()); // Normalize the normal

//             total_volume += (Q_F - CGAL::ORIGIN) * N_F * areaF;
//         }

//         return std::abs(total_volume) / 3.0;
//     }

// };


// double computeWeight(const Component& comp) {
//     int numTriangles = comp.triangles.size();
//     return 1.0 / (numTriangles * numTriangles);
// }

double Stitj(const Triangle& ti, const Triangle& tj, const std::vector<double>& weights) {
    double wA = weights[0], wL = weights[1], wS = weights[2], wN = weights[3];
    
    double areaDenominator = std::max(ti.area, tj.area);
    double maxLengthDenominator = ti.e_max + tj.e_max;
    double minLengthDenominator = ti.e_min + tj.e_min;

    double areaDifference = (areaDenominator != 0.0) ? wA * std::abs(ti.area - tj.area) / areaDenominator : 0.0;
    double maxLengthDifference = (maxLengthDenominator != 0.0) ? wL * std::abs(ti.e_max - tj.e_max) / maxLengthDenominator : 0.0;
    double minLengthDifference = (minLengthDenominator != 0.0) ? wS * std::abs(ti.e_min - tj.e_min) / minLengthDenominator : 0.0;
    double normalDifference = wN * (1 - (ti.normal.x * tj.normal.x + ti.normal.y * tj.normal.y + ti.normal.z * tj.normal.z)) / 2.0;

    return areaDifference + maxLengthDifference + minLengthDifference + normalDifference;
}



//MAIN THING COMES DOWN TO THIS DISJOINT SIMILAR TRIANGLES FUNCTION

std::vector<Triangle> findLargestDisjointSimilarTriangles(const std::vector<Triangle>& triangles) {
    std::vector<Triangle> seeds;

    // Sort triangles by area
    std::vector<Triangle> sortedTriangles = triangles;
    std::sort(sortedTriangles.begin(), sortedTriangles.end(), [](const Triangle& a, const Triangle& b) {
        return a.area > b.area;
    });

    // Set weights, with wN not equal to 0
    std::vector<double> weights = {.25, .25, .25 };
    double wN = 0.5;
    weights.push_back(wN);

    for (const Triangle& triangle : sortedTriangles) {
        bool isSimilar = false;
        for (const Triangle& seed : seeds) {
            if (Stitj(triangle, seed, weights) <= .01) { // Close to 0, meaning triangles are similar
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
            double similarity = Stitj(sortedTriangles[i], cluster.triangles[0], weights);
            if (similarity <= tau_S) {
                cluster.triangles.push_back(sortedTriangles[i]);
                assigned = true;
                break;
            }
        }
        if (!assigned) {
            // Create a new cluster for this triangle
            clusters.push_back({{sortedTriangles[i]}});
        }
    }
    return clusters;
}




double TriangleClusterSimilarity(const Triangle& triangle, const Cluster& cluster, const std::vector<double>& weights) {
    // Compute average properties for the cluster
    double avgArea = 0.0;
    Vector3D avgNormal;
    double avgE_min = 0.0;
    double avgE_max = 0.0;

    for (const Triangle& t : cluster.triangles) {
        avgArea += t.area;
        avgNormal += t.normal;
        avgE_min += t.e_min;
        avgE_max += t.e_max;
    }

    int clusterSize = cluster.triangles.size();
    avgArea /= clusterSize;
    avgNormal /= clusterSize;
    avgE_min /= clusterSize;
    avgE_max /= clusterSize;

    // Create an average triangle for the cluster
    Triangle avgTriangle;
    avgTriangle.area = avgArea;
    avgTriangle.normal = avgNormal;
    avgTriangle.e_min = avgE_min;
    avgTriangle.e_max = avgE_max;

    return Stitj(triangle, avgTriangle, weights);
}

bool isMoreSimilarToCluster(const Triangle& triangle, const Cluster& cluster1, const Cluster& cluster2, const std::vector<double>& weights) {
    double sim1 = TriangleClusterSimilarity(triangle, cluster1, weights);
    double sim2 = TriangleClusterSimilarity(triangle, cluster2, weights);

    return sim1 < sim2; // Note: Since STITJ is a difference, a smaller value indicates more similarity.
}


// Triangle chooseNextTriangleByHull(const std::vector<Component>& tempComponents, const std::vector<Triangle>& remainingTriangles) {
//     Triangle bestTriangle;
//     double bestScore = std::numeric_limits<double>::max();  // initialize with a high value

//     for (const Triangle& triangle : remainingTriangles) {
//         std::vector<double> volumeChanges;
//         for (const Component& component : tempComponents) {
//             Hull originalHull(component);  // Create Hull using the component's triangles
//             Hull tempHull = originalHull;  // Copy the original hull to a temporary one
            
//             tempHull.add(triangle);
//             tempHull.computeHull(); // Compute the convex hull before calculating its volume
//             volumeChanges.push_back(std::abs(tempHull.volume() - originalHull.volume()));
//         }

//         double meanVolumeChange = std::accumulate(volumeChanges.begin(), volumeChanges.end(), 0.0) / volumeChanges.size();
//         double variance = std::inner_product(volumeChanges.begin(), volumeChanges.end(), volumeChanges.begin(), 0.0) / volumeChanges.size() - meanVolumeChange * meanVolumeChange;

//         // Here you combine your volume change score and similarity measure.
//         // Adjust the weights as per your requirements.
//         double combinedScore = meanVolumeChange + variance;  // Just an example combination

//         if (combinedScore < bestScore) {
//             bestScore = combinedScore;
//             bestTriangle = triangle;
//         }
//     }

//     return bestTriangle;
// }


bool isSimilarTriangle(const Triangle& triangle, const Cluster& cluster, const std::vector<double>& weights, double tau_S) {
    for (const auto& t : cluster.triangles) {
        if (Stitj(triangle, t, weights) <= tau_S) {
            return true;
        }
    }
    return false;
}
enum class SeedingResult {
    PERFECT,
    UNDER_SEEDED,
    OVER_SEEDED
};

// SeedingResult determineSeeding(const std::vector<Triangle>& seeds, const std::vector<Component>& components) {
//     size_t seedsFound = 0;
//     for (const auto& seed : seeds) {
//         for (const auto& component : components) {
//             if (std::find(component.triangles.begin(), component.triangles.end(), seed) != component.triangles.end()) {
//                 seedsFound++;
//                 break;
//             }
//         }
//     }

//     if (seedsFound == seeds.size()) {
//         return SeedingResult::PERFECT;
//     } else if (seedsFound < seeds.size()) {
//         return SeedingResult::UNDER_SEEDED;
//     } else {
//         return SeedingResult::OVER_SEEDED;
//     }
// }

// bool areComponentsSimilar(const Component& a, const Component& b, double similarityThreshold) {
//     Hull hullA(a);
//     Hull hullB(b);
    
//     hullA.computeHull();
//     hullB.computeHull();
    
//     double volumeA = hullA.volume();
//     double volumeB = hullB.volume();

//     // Compare the volumes of the convex hulls of the two components
//     double relativeDifference = std::abs(volumeA - volumeB) / std::max(volumeA, volumeB);
//     return relativeDifference <= similarityThreshold;
// }

// void mergeSimilarComponents(std::vector<Component>& components, double similarityThreshold) {
//     for (size_t i = 0; i < components.size(); ++i) {
//         for (size_t j = i + 1; j < components.size();) {
//             if (areComponentsSimilar(components[i], components[j], similarityThreshold)) {  // This function needs to be implemented based on your criteria
//                 components[i].triangles.insert(components[i].triangles.end(),
//                                                components[j].triangles.begin(),
//                                                components[j].triangles.end());
//                 components.erase(components.begin() + j);
//             } else {
//                 ++j;
//             }
//         }
//     }
// }

// bool recursiveBacktracking(const std::vector<Triangle>& remainingTriangles, std::vector<Component>& solution, const std::vector<double>& weights, double tau_S, int depth = 0, const int MAX_DEPTH = 3) {
//     if (depth > MAX_DEPTH) {
//         return false;  // Maximum depth reached without a solution.
//     }

//     if (remainingTriangles.empty()) {
//         return true;  // Found a solution.
//     }

//     for (const Triangle& triangle : remainingTriangles) {
//         for (Component& component : solution) {
//             Cluster tempCluster = {component.triangles};
//             if (isSimilarTriangle(triangle, tempCluster, weights, tau_S)) {
//                 component.triangles.push_back(triangle);
//                 std::vector<Triangle> newRemaining = remainingTriangles;
//                 newRemaining.erase(std::remove(newRemaining.begin(), newRemaining.end(), triangle), newRemaining.end());
//                 if (recursiveBacktracking(newRemaining, solution, weights, tau_S, depth + 1)) {
//                     return true;
//                 }
//                 component.triangles.pop_back();
//             }
//         }
//     }

//     return false;
// }


// std::vector<Component> randomizedGrowthOptimization(const Cluster& cluster) {
//     std::vector<Component> components;
//     std::vector<Triangle> remainingTriangles = cluster.triangles;
    
//     int iteration = 1; // To keep track of iterations
    
//     while (!remainingTriangles.empty()) {
//         std::cout << "Iteration " << iteration++ << std::endl;
        
//         // Seed with the largest disjoint similar triangles
//         std::vector<Triangle> seeds = findLargestDisjointSimilarTriangles(remainingTriangles);
//         std::cout << "Number of seeds found: " << seeds.size() << std::endl;
        
//         SeedingResult seedingResult = determineSeeding(seeds, components);
//         std::vector<Component> tempComponents(seeds.size());
//         for (size_t i = 0; i < seeds.size(); i++) {
//             tempComponents[i].triangles.push_back(seeds[i]);
//             remainingTriangles.erase(std::remove(remainingTriangles.begin(), remainingTriangles.end(), seeds[i]), remainingTriangles.end());
//         }

//         bool canGrow = true;
//         int backtrackSteps = 0;
//         const int MAX_BACKTRACK = 3;

//         while (canGrow && backtrackSteps < MAX_BACKTRACK) {
//             canGrow = false;
//             int grownTriangles = 0;
//             for (auto& component : tempComponents) {
//                 Triangle nextTriangle = chooseNextTriangleByHull(tempComponents, remainingTriangles);
//                 if (nextTriangle.isValid()) {
//                     component.triangles.push_back(nextTriangle);
//                     remainingTriangles.erase(std::remove(remainingTriangles.begin(), remainingTriangles.end(), nextTriangle), remainingTriangles.end());
//                     canGrow = true;
//                     grownTriangles++;
//                 } else {
//                     // If unable to grow, we might need to backtrack
//                     if (!component.triangles.empty()) {
//                         remainingTriangles.push_back(component.triangles.back());
//                         component.triangles.pop_back();
//                         backtrackSteps++;
//                     }
//                 }
//             }
//             std::cout << "Grown triangles in this step: " << grownTriangles << std::endl;
//         }
//         std::cout << "Total backtracks in this iteration: " << backtrackSteps << std::endl;

//         // Merge valid grown components
//         int validComponentsCount = 0;
//         for (const auto& comp : tempComponents) {
//             if (comp.isValid()) {
//                 components.push_back(comp);
//                 validComponentsCount++;
//             }
//         }
//         std::cout << "Valid components in this iteration: " << validComponentsCount << std::endl;
//     }

//     double similarityThreshold = 0.9;  // or your desired threshold
//     mergeSimilarComponents(components, similarityThreshold);

//     std::cout << "Total components after merging: " << components.size() << std::endl;

//     return components;
// }


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

void refineClusters(std::vector<Cluster>& clusters, double tau_N) {
    int iteration = 0;
    SpatialHash spatialHash(0.1);
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

    size_t initialTriangleCount = 0;
    for (const auto& cluster : clusters) {
        initialTriangleCount += cluster.triangles.size();
    }

    while (iteration < 10) {
        std::cout << "Starting iteration " << iteration << " with " << clusters.size() << " clusters." << std::endl;

        int mergesInThisIteration = 0;
        int splitsInThisIteration = 0;
        size_t pairsProcessed = 0;  // Initialize counter for this iteration

        // Declare a variable to accumulate the total time for the specific loop
        std::chrono::milliseconds totalLoopTime(0);

        std::map<std::pair<size_t, size_t>, std::vector<Triangle>> mergeMap;
        std::map<std::pair<size_t, size_t>, std::vector<Triangle>> splitMap;
        std::unordered_map<size_t, std::pair<size_t, double>> mostSimilarCluster; // sourceCluster -> (targetCluster, similarity)

        for (size_t i = 0; i < clusters.size(); ++i) {
            std::unordered_set<Triangle> uniqueTriangles(clusters[i].triangles.begin(), clusters[i].triangles.end());
            if (uniqueTriangles.size() != clusters[i].triangles.size()) {
                std::cout << "Duplicate triangles found in cluster " << i << std::endl;
            }

            if (clusters[i].triangles.empty()) continue;

            std::unordered_set<int> neighboringClusterIndices = spatialHash.getNeighboringClustersForCluster(clusters[i]);
            std::unordered_map<Triangle, std::vector<Triangle>> potentialNeighborsCache;

            for (const auto& triangleA : clusters[i].triangles) {
                potentialNeighborsCache[triangleA] = spatialHash.getPotentialNeighbors(triangleA);
            }

            std::vector<int> neighboringClusterIndicesVec(neighboringClusterIndices.begin(), neighboringClusterIndices.end());

            auto loopStartTime = std::chrono::high_resolution_clock::now();  // Start time measurement
            std::unordered_set<Triangle> alreadyAddedToMaps;  // New line


            #pragma omp parallel for  // Enable OpenMP parallelization
            for (size_t idx = 0; idx < neighboringClusterIndicesVec.size(); ++idx) {
                int j = neighboringClusterIndicesVec[idx];
                if (i == j || clusters[j].triangles.empty()) continue;

                std::unordered_set<Triangle> localMergeList;
                std::vector<Triangle> localSplitListI;
                std::vector<Triangle> localSplitListJ;
                std::unordered_set<Triangle> tbTriangles(clusters[j].triangles.begin(), clusters[j].triangles.end());

                double similarity = ATaTb(clusters[i], clusters[j], potentialNeighborsCache, triangleAdjacency, tbTriangles);
                std::unordered_set<Triangle> alreadyProcessed;
                std::unordered_set<Triangle> toBeMerged;
                std::unordered_set<Triangle> toBeSplit;

                if (similarity >= tau_N) {
                    for (const auto& triangleA : clusters[i].triangles) {
                        if (alreadyProcessed.count(triangleA) > 0) continue;  // Skip if already processed

                        bool foundAdjacent = false;
                        for (const auto& triangleB : potentialNeighborsCache[triangleA]) {
                            if (triangleA.isAdjacent(triangleB, triangleAdjacency)) {
                                // Both triangles meet the merge condition, so add them to the list
                                toBeMerged.insert(triangleA);
                                toBeMerged.insert(triangleB);  // Insert triangleB as well
                                alreadyProcessed.insert(triangleA);  // Mark as processed
                                alreadyProcessed.insert(triangleB);  // Mark as processed
                                foundAdjacent = true;
                                break;
                            }
                        }
                        if (!foundAdjacent) {
                            toBeSplit.insert(triangleA);
                            alreadyProcessed.insert(triangleA);  // Mark as processed
                        }
                    }

                    for (const auto& triangleB : clusters[j].triangles) {
                        if (alreadyProcessed.count(triangleB) > 0) continue;  // Skip if already processed
                        if (toBeMerged.count(triangleB) == 0) {
                            // Only add to be split if it's not in the merge list
                            toBeSplit.insert(triangleB);
                            alreadyProcessed.insert(triangleB);  // Mark as processed
                        }
                    }

                    #pragma omp critical  // Critical section
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

        // Create a temporary copy of the clusters
        std::vector<Cluster> tempClusters = clusters;

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
        }

        clusters = tempClusters;  // Update clusters with merged state

        std::unordered_set<Triangle> splitTriangles;  // To keep track of triangles that have been split
        std::chrono::microseconds totalTime(0);
        std :: cout << "cluster size before split: " << clusters.size() << std::endl;

        tempClusters = clusters;  // Start with the most current data

        // Loop through each entry in the splitMap
        for (auto& entry : splitMap) {
            auto& keyPair = entry.first;
            auto& trianglesToSplit = entry.second;
            size_t sourceClusterIndex = keyPair.first;
            size_t targetClusterIndex = keyPair.second;

            // Insert initial count checks here
            size_t initialSourceCount = tempClusters[sourceClusterIndex].triangles.size();
            size_t initialTargetCount = tempClusters[targetClusterIndex].triangles.size();
            size_t initialTotalCount = initialSourceCount + initialTargetCount;


            auto splitStart = std::chrono::high_resolution_clock::now();  // Start timing

            std::vector<std::unordered_set<Triangle>> privateSplitTriangles(omp_get_max_threads());
            std::vector<std::vector<Triangle>> privateSourceTriangles(omp_get_max_threads());
            std::vector<std::vector<Triangle>> privateTargetTriangles(omp_get_max_threads());

            // Parallel loop for splitting
            #pragma omp parallel for
            for (size_t i = 0; i < trianglesToSplit.size(); ++i) {
                const auto& triangle = trianglesToSplit[i];
                int thread_id = omp_get_thread_num();
                
                // Skip if already split
                if (privateSplitTriangles[thread_id].find(triangle) != privateSplitTriangles[thread_id].end()) {
                    continue;
                }

                auto& sourceTriangles = tempClusters[sourceClusterIndex].triangles;
                auto& targetTriangles = tempClusters[targetClusterIndex].triangles;
                
                if (std::find(sourceTriangles.begin(), sourceTriangles.end(), triangle) != sourceTriangles.end()) {
                    privateSourceTriangles[thread_id].push_back(triangle);
                } else if (std::find(targetTriangles.begin(), targetTriangles.end(), triangle) != targetTriangles.end()) {
                    privateTargetTriangles[thread_id].push_back(triangle);
                }

                privateSplitTriangles[thread_id].insert(triangle);
            }
            Cluster newClusterSource, newClusterTarget;
            // Merging thread-private data back to main structures
            std::vector<Triangle> mergedSourceTriangles, mergedTargetTriangles;
            for (int i = 0; i < omp_get_max_threads(); ++i) {
                mergedSourceTriangles.insert(mergedSourceTriangles.end(), privateSourceTriangles[i].begin(), privateSourceTriangles[i].end());
                mergedTargetTriangles.insert(mergedTargetTriangles.end(), privateTargetTriangles[i].begin(), privateTargetTriangles[i].end());
            }

            // Now you can safely remove these triangles from the main clusters
            // This part is not parallelized as it involves erasing elements
            for (const auto& triangle : mergedSourceTriangles) {
                auto itSource = std::remove(tempClusters[sourceClusterIndex].triangles.begin(), tempClusters[sourceClusterIndex].triangles.end(), triangle);
                tempClusters[sourceClusterIndex].triangles.erase(itSource, tempClusters[sourceClusterIndex].triangles.end());
            }

            for (const auto& triangle : mergedTargetTriangles) {
                auto itTarget = std::remove(tempClusters[targetClusterIndex].triangles.begin(), tempClusters[targetClusterIndex].triangles.end(), triangle);
                tempClusters[targetClusterIndex].triangles.erase(itTarget, tempClusters[targetClusterIndex].triangles.end());
            }

            auto splitEnd = std::chrono::high_resolution_clock::now();  // Stop timing
            auto splitDuration = std::chrono::duration_cast<std::chrono::microseconds>(splitEnd - splitStart);
            totalTime += splitDuration;

            // Populate the new clusters with the merged triangles
            newClusterSource.triangles = mergedSourceTriangles;
            newClusterTarget.triangles = mergedTargetTriangles;

            // Now you can safely remove these triangles from the main clusters
            // This part is not parallelized as it involves erasing elements
            for (const auto& triangle : mergedSourceTriangles) {
                auto itSource = std::remove(tempClusters[sourceClusterIndex].triangles.begin(), tempClusters[sourceClusterIndex].triangles.end(), triangle);
                tempClusters[sourceClusterIndex].triangles.erase(itSource, tempClusters[sourceClusterIndex].triangles.end());
            }

            for (const auto& triangle : mergedTargetTriangles) {
                auto itTarget = std::remove(tempClusters[targetClusterIndex].triangles.begin(), tempClusters[targetClusterIndex].triangles.end(), triangle);
                tempClusters[targetClusterIndex].triangles.erase(itTarget, tempClusters[targetClusterIndex].triangles.end());
            }

            // Insert final count checks and validation here
            size_t finalSourceCount = tempClusters[sourceClusterIndex].triangles.size();
            size_t finalTargetCount = tempClusters[targetClusterIndex].triangles.size();
            size_t finalTotalCount = finalSourceCount + finalTargetCount;

            // Add the new clusters
            if (!newClusterSource.triangles.empty()) {
                tempClusters.push_back(newClusterSource);
            }
            if (!newClusterTarget.triangles.empty()) {
                tempClusters.push_back(newClusterTarget);
            }
        }
        
        clusters = tempClusters;

        validateTotalTriangles(tempClusters, initialTriangleCount);

        std::cout << "Total time taken: " << totalTime.count() << " microseconds" << std::endl;

        std::cout << "After " << iteration << " iterations, merged " << mergeMap.size() << " triangles and split " << splitMap.size() << " triangles." << std::endl;
        // clusters.insert(clusters.end(), newClusters.begin(), newClusters.end());


        //clusters before erase
        std::cout << "Before erase: " << clusters.size() << std::endl;

        // Remove empty clusters after merges and splits
        auto it = std::remove_if(clusters.begin(), clusters.end(),
            [](const Cluster& cluster) { return cluster.triangles.empty(); });
            
        std::cout << "Number of empty clusters: " << std::distance(it, clusters.end()) << std::endl;

        clusters.erase(it, clusters.end());
        
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

        // for (const auto& [triangle, count] : triangleCount) {
        //     if (count > 1) {
        //         std::cout << "Triangle " << triangle.toString() << " appears " << count << " times." << std::endl;
        //     }
        // }



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


std::vector<Cluster> createSearchSpaces(const std::vector<Triangle>& triangles, const std::vector<double>& weights, double tau_S) {
    // 1. Initial clustering based on shape similarity
    std::cout << "1. Initial clustering based on shape similarity" << std::endl;
    std::vector<Cluster> clusters = initialClusteringByShapeSimilarity(triangles, weights, tau_S);
    std::cout << "Number of initial clusters: " << clusters.size() << std::endl;
    std::cout << "Initial clustering done. Number of clusters: " << clusters.size() << std::endl;
    
    // 3. Iterative merging and splitting based on adjacency similarity
    std::cout << "3. Iterative merging and splitting based on adjacency similarity" << std::endl;
    double tau_N = 0.05;  // Threshold for adjacency similarity
    refineClusters(clusters, tau_N);

    // std::set<Triangle> usedTriangles; // This will store triangles that are part of a component
    // std::vector<Component> allComponents; // This will store all grown components

    // // 4. Randomized Growth Optimization
    // std::cout << "4. Randomized Growth Optimization" << std::endl;
    // for (Cluster& cluster : clusters) {
    //     std::vector<Component> components = randomizedGrowthOptimization(cluster);
    //     allComponents.insert(allComponents.end(), components.begin(), components.end());

    //     // For each component, add its triangles to the usedTriangles set
    //     for (const Component& comp : components) {
    //         usedTriangles.insert(comp.triangles.begin(), comp.triangles.end());
    //     }
    // }

    // std::vector<Triangle> remainingTriangles;
    // for (const Triangle& tri : triangles) {
    //     if (usedTriangles.find(tri) == usedTriangles.end()) {
    //         remainingTriangles.push_back(tri);
    //     }
    // }

    // // 5. Recursive Backtracking
    // std::cout << "5. Recursive Backtracking" << std::endl;
    // bool success = recursiveBacktracking(remainingTriangles, allComponents, weights, tau_S);
    // if (!success) {
    //     std::cout << "Some triangles couldn't be assigned to any component! Re-seeding..." << std::endl;

    //     // Re-seed with leftover triangles and start another iteration of randomized growth
    //     for (Cluster& cluster : clusters) {
    //         std::vector<Triangle> seeds = findLargestDisjointSimilarTriangles(cluster.triangles);
    //         cluster.triangles = seeds;  // Replace the cluster's triangles with the seeds
    //         randomizedGrowthOptimization(cluster);
    //     }

    //     // Attempt recursive backtracking again
    //     success = recursiveBacktracking(remainingTriangles, allComponents, weights, tau_S);
    // }


    // Save initial clusters to an OBJ file for visualization
    float hueStep = 360.0f / clusters.size();
    float currentHue = 0.0f;
    float saturation = 1.0f;
    float value = 1.0f;

    for (auto& cluster : clusters) {
        Vector3D assignedColor = HSVtoRGB(currentHue, saturation, value);
        for (auto& triangle : cluster.triangles) {
            triangle.setColor(assignedColor);
        }
        currentHue += hueStep;
    }

    std::cout << "Number of clusters after refining: " << clusters.size() << std::endl;

    // Save the clusters after refining
    saveToOBJ(clusters, "clusters");
    std::cout << "Saved refined clusters to OBJ for visualization." << std::endl;

    return clusters;
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




int main() {
    std::cout << "Starting program..." << std::endl;
    // std::vector<Triangle> allTriangles = {};  // Fill with triangles from your data source
    std::vector<double> weights = {1.0, 1.0, 1.0, 1.0};
    std::vector<Vector3D> allVertices;    
    std::vector<std::vector<Vector3D>> quadPolygons; // Store polygons with 4 vertices
    double tau_S = 0.01;  // Threshold for shape similarity
    std::string inputfile = "/home/kingpin/Documents/blender/kb3d_americana-native.obj";
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