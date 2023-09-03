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
class Triangle; // forward declare the Triangle class
bool areAdjacent(const Triangle& t1, const Triangle& t2);
class Component; // forward declare the Component class
bool areTrianglesContiguous(const std::vector<Triangle>& triangles);
bool areComponentsSimilar(const Component& a, const Component& b, double similarityThreshold);
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

class Triangle {
public:
    Vector3D vec1, vec2, vec3;  // Geometric vertices
    Vector3D normal;
    double area;
    double e_min;
    double e_max;
    Vector3D color;  // RGB color for the triangle

    // Default constructor
    Triangle() : area(0), e_min(0), e_max(0) {}

    // Constructor using Vector3D
    Triangle(const Vector3D& a, const Vector3D& b, const Vector3D& c)
        : vec1(a), vec2(b), vec3(c) {
        // Compute edges
        Vector3D edge1 = vec2 - vec1;
        Vector3D edge2 = vec3 - vec1;

        // Compute normal
        normal = edge1.cross(edge2);
        double normalLength = normal.length();
        normal = normal / normalLength;  // Normalize the normal vector

        // Compute area
        area = 0.5 * normalLength;

        // Compute edge lengths
        double len1 = edge1.length();
        double len2 = edge2.length();
        double len3 = (vec3 - vec2).length();
        e_min = std::min({len1, len2, len3});
        e_max = std::max({len1, len2, len3});
    }

    // Check if triangle is valid based on vertices
    bool isValid() const {
        return !vec1.isZero() || !vec2.isZero() || !vec3.isZero();
    }

    bool operator==(const Triangle& other) const {
        return (vec1 == other.vec1 && vec2 == other.vec2 && vec3 == other.vec3) ||
               (vec1 == other.vec1 && vec2 == other.vec3 && vec3 == other.vec2) ||
               (vec1 == other.vec2 && vec2 == other.vec1 && vec3 == other.vec3) ||
               (vec1 == other.vec2 && vec2 == other.vec3 && vec3 == other.vec1) ||
               (vec1 == other.vec3 && vec2 == other.vec1 && vec3 == other.vec2) ||
               (vec1 == other.vec3 && vec2 == other.vec2 && vec3 == other.vec1);
    }

    bool operator!=(const Triangle& other) const {
        return !(*this == other);
    }

    bool operator<(const Triangle& other) const {
        return this->area < other.area;
    }

    bool isAdjacent(const Triangle& other) const {
        return areAdjacent(*this, other);
    }

    void setColor(const Vector3D& assignedColor) {
        color = assignedColor;
    }

    std::string toString() const {
        return "[ " + vec1.toString() + ", " + vec2.toString() + ", " + vec3.toString() + " ]";
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

        // Print the vertices to verify the mapping
        // std::cout << "Mapping Triangle Vertices: (" << v1.x << "," << v1.y << "," << v1.z << "), "
        //           << "(" << v2.x << "," << v2.y << "," << v2.z << "), "
        //           << "(" << v3.x << "," << v3.y << "," << v3.z << ")" << std::endl;

        Triangle tri(v1, v2, v3);
        triangles.push_back(tri);
    }
    return triangles;
}

struct Edge {
    int v1, v2;

    // Constructor to ensure the smaller vertex is always first. This helps in checking for duplicate edges.
    Edge(int a, int b) {
        if (a < b) {
            v1 = a;
            v2 = b;
        } else {
            v1 = b;
            v2 = a;
        }
    }

    // We need to provide a custom comparator for the map
    bool operator<(const Edge& other) const {
        if (v1 == other.v1) return v2 < other.v2;
        return v1 < other.v1;
    }
};


struct EdgeIndices {
    int v1, v2;

    bool operator==(const EdgeIndices& other) const {
        return (v1 == other.v1 && v2 == other.v2) || (v1 == other.v2 && v2 == other.v1);
    }
};

namespace std {
    template<>
    struct hash<EdgeIndices> {
        size_t operator()(const EdgeIndices& edge) const {
            return std::hash<int>()(edge.v1) ^ std::hash<int>()(edge.v2);
        }
    };
}

struct TriangleHash {
    std::size_t operator()(const Triangle& t) const {
        std::hash<Vertex> vertexHasher;
        size_t h1 = vertexHasher(t.vec1);
        size_t h2 = vertexHasher(t.vec2);
        size_t h3 = vertexHasher(t.vec3);
        return h1 ^ (h2 << 1) ^ (h3 << 2);  // or use boost::hash_combine (better)
        }
};

namespace std {
    template <>
    struct hash<Triangle> {
        std::size_t operator()(const Triangle& t) const {
            TriangleHash triangleHasher;
            return triangleHasher(t);
        }
    };
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
        return std::get<0>(k) ^ std::get<1>(k) ^ std::get<2>(k);
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
        for (const Vector3D& vec : {t.vec1, t.vec2, t.vec3}) {
            hashTable[hash(vec)].emplace_back(IndexedTriangle{clusterIndex, t});
        }
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
        
        // Cache hash values
        std::array<std::tuple<int, int, int>, 3> hashes = {hash(t.vec1), hash(t.vec2), hash(t.vec3)};
        
        for (const auto& h : hashes) {
            for (const auto& indexedTriangle : hashTable[h]) {
                clusterIndices.insert(indexedTriangle.clusterIndex);
            }
        }
        
        return clusterIndices;
    }
    std::vector<Triangle> getPotentialNeighbors(const Triangle& t) {
        std::vector<Triangle> neighbors;

        // Cache the hashes of vertices
        std::array<std::tuple<int, int, int>, 3> hashes = {hash(t.vec1), hash(t.vec2), hash(t.vec3)};
            
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


    //original
    std::unordered_set<int> getNeighboringClustersForCluster(const Cluster& cluster) {
        std::unordered_set<int> clusterIndices;
        std::unordered_set<Triangle> processedTriangles;
            
        for (const auto& triangle : cluster.triangles) {
            if (!processedTriangles.insert(triangle).second) continue;
            auto neighboringClustersForTriangle = getNeighboringClusters(triangle);
            clusterIndices.insert(neighboringClustersForTriangle.begin(), neighboringClustersForTriangle.end());
        }
            
        return clusterIndices;
    }

    //original
    std::vector<Triangle> getPotentialNeighborsForCluster(const Cluster& cluster) {
        std::vector<Triangle> neighbors;
        std::unordered_set<Triangle> processedTriangles;
            
        for (const auto& triangle : cluster.triangles) {
            if (!processedTriangles.insert(triangle).second) continue;
            auto potentialNeighbors = getPotentialNeighbors(triangle);
            neighbors.insert(neighbors.end(), potentialNeighbors.begin(), potentialNeighbors.end());
        }

        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());

        return neighbors;
    }



    void clear() {
        hashTable.clear();
    }
};

void preprocessSpatialHash(const std::vector<Cluster>& clusters, SpatialHash& spatialHash) {
    for (size_t clusterIndex = 0; clusterIndex < clusters.size(); ++clusterIndex) {
        const Cluster& cluster = clusters[clusterIndex];
        for (const Triangle& triangle : cluster.triangles) {
            spatialHash.insert(triangle, clusterIndex);
        }
    }
}

std::vector<Triangle> getAdjacentTriangles(const Triangle& target, const std::vector<Triangle>& triangles) {
    std::vector<Triangle> adjacentTriangles;

    for (const Triangle& triangle : triangles) {
        if (areAdjacent(target, triangle)) {  // This function checks if two triangles share an edge or vertex
            adjacentTriangles.push_back(triangle);
        }
    }

    return adjacentTriangles;
}

bool areTrianglesContiguous(const std::vector<Triangle>& triangles) {
    if (triangles.empty()) return false;

    std::unordered_set<Triangle> visited;
    std::stack<Triangle> stack;

    stack.push(triangles[0]);

    while (!stack.empty()) {
        Triangle current = stack.top();
        stack.pop();

        if (visited.find(current) == visited.end()) {
            visited.insert(current);

            for (Triangle adjacent : getAdjacentTriangles(current, triangles)) {  
                if (visited.find(adjacent) == visited.end()) {
                    stack.push(adjacent);
                }
            }
        }
    }

    return visited.size() == triangles.size();
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



void saveTrianglesOnlyToOBJ(const std::vector<Triangle>& triangles, const std::string& filename) {
    std::ostringstream objStream;

    // Generate vertex and face data
    for (const auto& triangle : triangles) {
        objStream << "v " << triangle.vec1.x << " " << triangle.vec1.y << " " << triangle.vec1.z << "\n";
        objStream << "v " << triangle.vec2.x << " " << triangle.vec2.y << " " << triangle.vec2.z << "\n";
        objStream << "v " << triangle.vec3.x << " " << triangle.vec3.y << " " << triangle.vec3.z << "\n";
    }
    for (size_t i = 0; i < triangles.size(); ++i) {
        size_t baseIdx = 1 + i * 3;
        objStream << "f " << baseIdx << " " << (baseIdx + 1) << " " << (baseIdx + 2) << "\n";
    }

    // Check if there's an issue with the stream
    if (objStream.fail()) {
        std::cerr << "Error: Failed to write to the stream." << std::endl;
        return;
    }

    // Open the .obj file and check if it opened correctly
    std::ofstream objFile(filename + ".obj");
    if (!objFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing: " << (filename + ".obj") << std::endl;
        return;
    }

    // Write to the file, flush the stream, and then close it
    objFile << objStream.str();
    objFile.flush();
    objFile.close();

    // Print debug information
    std::cout << "Saved " << triangles.size() << " triangles to " << (filename + ".obj") << std::endl;
}



// bool areAdjacent(const Triangle& t1, const Triangle& t2) {
//     // Helper function to check if a triangle contains a vertex (Vector3D in this case)
//     auto triangleContains = [](const Triangle& t, const Vector3D& v) {
//         return t.vec1 == v || t.vec2 == v || t.vec3 == v;
//     };

//     if (t1 == t2) {
//         return false;
//     }

//     // Check if the triangles share any vertices
//     int sharedVerticesCount = 0;
//     if (triangleContains(t2, t1.vec1)) {
//         sharedVerticesCount++;
//     }
//     if (triangleContains(t2, t1.vec2)) {
//         sharedVerticesCount++;
//     }
//     if (triangleContains(t2, t1.vec3)) {
//         sharedVerticesCount++;
//     }

//     // If they share exactly two vertices, they are adjacent
//     if (sharedVerticesCount == 2) {
//         return true;
//     } else {
//         return false;
//     }
// }

bool areAdjacent(const Triangle& t1, const Triangle& t2) {
    if (t1 == t2) {
        return false;
    }

    // Put vertices of t2 in a set for faster lookup
    std::unordered_set<Vector3D, Vector3DHash> t2Vertices = {t2.vec1, t2.vec2, t2.vec3};

    // Count the vertices of t1 that are in t2
    int sharedVerticesCount = 0;
    sharedVerticesCount += t2Vertices.count(t1.vec1);
    sharedVerticesCount += t2Vertices.count(t1.vec2);
    sharedVerticesCount += t2Vertices.count(t1.vec3);

    // If they share exactly two vertices, they are adjacent
    return sharedVerticesCount == 2;
}




struct Component {
    std::vector<Triangle> triangles;

    bool isValid() const {
        // Check if triangles are contiguous
        // This requires a function that checks if a group of triangles form a connected set
        if (!areTrianglesContiguous(triangles)) {  // To be implemented
            return false;
        }
        // Check other properties (iii and iv) if needed
        return true;
    }
    
    double weight() const {
        return 1.0 / (triangles.size() * triangles.size());
    }
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
    }


    // Add a triangle to the collection of points
    void add(const Triangle& triangle) {
        points.push_back(Point_3(triangle.vec1.x, triangle.vec1.y, triangle.vec1.z));
        points.push_back(Point_3(triangle.vec2.x, triangle.vec2.y, triangle.vec2.z));
        points.push_back(Point_3(triangle.vec3.x, triangle.vec3.y, triangle.vec3.z));
    }

    // Compute the convex hull and store it in Polyhedron P
    void computeHull() {
        CGAL::convex_hull_3(points.begin(), points.end(), P);
    }

    double volume() const {
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

        return std::abs(total_volume) / 3.0;
    }

};


double computeWeight(const Component& comp) {
    int numTriangles = comp.triangles.size();
    return 1.0 / (numTriangles * numTriangles);
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

    return areaDifference + maxLengthDifference + minLengthDifference + normalDifference;
}

double ATaTb(const Cluster& Ta, const Cluster& Tb, SpatialHash& spatialHash) {
    int count = 0;

    // Convert Tb's triangles to a set for O(1) lookup
    std::unordered_set<Triangle> tbTriangles(Tb.triangles.begin(), Tb.triangles.end());

    // For each triangle in Ta, check its potential neighbors
    for (const auto& ta_i : Ta.triangles) {
        auto potentialNeighbors = spatialHash.getPotentialNeighbors(ta_i);
        for (const auto& potentialNeighbor : potentialNeighbors) {
            // Only check adjacency if the potential neighbor is in Tb
            if (tbTriangles.count(potentialNeighbor) && ta_i.isAdjacent(potentialNeighbor)) {
                count++;
            }
        }
    }

    double result = static_cast<double>(count) / std::min(Ta.triangles.size(), Tb.triangles.size());
    return result;
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


Triangle chooseNextTriangleByHull(const std::vector<Component>& tempComponents, const std::vector<Triangle>& remainingTriangles) {
    Triangle bestTriangle;
    double bestScore = std::numeric_limits<double>::max();  // initialize with a high value

    for (const Triangle& triangle : remainingTriangles) {
        std::vector<double> volumeChanges;
        for (const Component& component : tempComponents) {
            Hull originalHull(component);  // Create Hull using the component's triangles
            Hull tempHull = originalHull;  // Copy the original hull to a temporary one
            
            tempHull.add(triangle);
            tempHull.computeHull(); // Compute the convex hull before calculating its volume
            volumeChanges.push_back(std::abs(tempHull.volume() - originalHull.volume()));
        }

        double meanVolumeChange = std::accumulate(volumeChanges.begin(), volumeChanges.end(), 0.0) / volumeChanges.size();
        double variance = std::inner_product(volumeChanges.begin(), volumeChanges.end(), volumeChanges.begin(), 0.0) / volumeChanges.size() - meanVolumeChange * meanVolumeChange;

        // Here you combine your volume change score and similarity measure.
        // Adjust the weights as per your requirements.
        double combinedScore = meanVolumeChange + variance;  // Just an example combination

        if (combinedScore < bestScore) {
            bestScore = combinedScore;
            bestTriangle = triangle;
        }
    }

    return bestTriangle;
}


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

SeedingResult determineSeeding(const std::vector<Triangle>& seeds, const std::vector<Component>& components) {
    size_t seedsFound = 0;
    for (const auto& seed : seeds) {
        for (const auto& component : components) {
            if (std::find(component.triangles.begin(), component.triangles.end(), seed) != component.triangles.end()) {
                seedsFound++;
                break;
            }
        }
    }

    if (seedsFound == seeds.size()) {
        return SeedingResult::PERFECT;
    } else if (seedsFound < seeds.size()) {
        return SeedingResult::UNDER_SEEDED;
    } else {
        return SeedingResult::OVER_SEEDED;
    }
}

bool areComponentsSimilar(const Component& a, const Component& b, double similarityThreshold) {
    Hull hullA(a);
    Hull hullB(b);
    
    hullA.computeHull();
    hullB.computeHull();
    
    double volumeA = hullA.volume();
    double volumeB = hullB.volume();

    // Compare the volumes of the convex hulls of the two components
    double relativeDifference = std::abs(volumeA - volumeB) / std::max(volumeA, volumeB);
    return relativeDifference <= similarityThreshold;
}

void mergeSimilarComponents(std::vector<Component>& components, double similarityThreshold) {
    for (size_t i = 0; i < components.size(); ++i) {
        for (size_t j = i + 1; j < components.size();) {
            if (areComponentsSimilar(components[i], components[j], similarityThreshold)) {  // This function needs to be implemented based on your criteria
                components[i].triangles.insert(components[i].triangles.end(),
                                               components[j].triangles.begin(),
                                               components[j].triangles.end());
                components.erase(components.begin() + j);
            } else {
                ++j;
            }
        }
    }
}

bool recursiveBacktracking(const std::vector<Triangle>& remainingTriangles, std::vector<Component>& solution, const std::vector<double>& weights, double tau_S, int depth = 0, const int MAX_DEPTH = 3) {
    if (depth > MAX_DEPTH) {
        return false;  // Maximum depth reached without a solution.
    }

    if (remainingTriangles.empty()) {
        return true;  // Found a solution.
    }

    for (const Triangle& triangle : remainingTriangles) {
        for (Component& component : solution) {
            Cluster tempCluster = {component.triangles};
            if (isSimilarTriangle(triangle, tempCluster, weights, tau_S)) {
                component.triangles.push_back(triangle);
                std::vector<Triangle> newRemaining = remainingTriangles;
                newRemaining.erase(std::remove(newRemaining.begin(), newRemaining.end(), triangle), newRemaining.end());
                if (recursiveBacktracking(newRemaining, solution, weights, tau_S, depth + 1)) {
                    return true;
                }
                component.triangles.pop_back();
            }
        }
    }

    return false;
}


std::vector<Component> randomizedGrowthOptimization(const Cluster& cluster) {
    std::vector<Component> components;
    std::vector<Triangle> remainingTriangles = cluster.triangles;
    
    int iteration = 1; // To keep track of iterations
    
    while (!remainingTriangles.empty()) {
        std::cout << "Iteration " << iteration++ << std::endl;
        
        // Seed with the largest disjoint similar triangles
        std::vector<Triangle> seeds = findLargestDisjointSimilarTriangles(remainingTriangles);
        std::cout << "Number of seeds found: " << seeds.size() << std::endl;
        
        SeedingResult seedingResult = determineSeeding(seeds, components);
        std::vector<Component> tempComponents(seeds.size());
        for (size_t i = 0; i < seeds.size(); i++) {
            tempComponents[i].triangles.push_back(seeds[i]);
            remainingTriangles.erase(std::remove(remainingTriangles.begin(), remainingTriangles.end(), seeds[i]), remainingTriangles.end());
        }

        bool canGrow = true;
        int backtrackSteps = 0;
        const int MAX_BACKTRACK = 3;

        while (canGrow && backtrackSteps < MAX_BACKTRACK) {
            canGrow = false;
            int grownTriangles = 0;
            for (auto& component : tempComponents) {
                Triangle nextTriangle = chooseNextTriangleByHull(tempComponents, remainingTriangles);
                if (nextTriangle.isValid()) {
                    component.triangles.push_back(nextTriangle);
                    remainingTriangles.erase(std::remove(remainingTriangles.begin(), remainingTriangles.end(), nextTriangle), remainingTriangles.end());
                    canGrow = true;
                    grownTriangles++;
                } else {
                    // If unable to grow, we might need to backtrack
                    if (!component.triangles.empty()) {
                        remainingTriangles.push_back(component.triangles.back());
                        component.triangles.pop_back();
                        backtrackSteps++;
                    }
                }
            }
            std::cout << "Grown triangles in this step: " << grownTriangles << std::endl;
        }
        std::cout << "Total backtracks in this iteration: " << backtrackSteps << std::endl;

        // Merge valid grown components
        int validComponentsCount = 0;
        for (const auto& comp : tempComponents) {
            if (comp.isValid()) {
                components.push_back(comp);
                validComponentsCount++;
            }
        }
        std::cout << "Valid components in this iteration: " << validComponentsCount << std::endl;
    }

    double similarityThreshold = 0.9;  // or your desired threshold
    mergeSimilarComponents(components, similarityThreshold);

    std::cout << "Total components after merging: " << components.size() << std::endl;

    return components;
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

int countAdjacentTriangles(const Cluster& Ta, const Cluster& Tb) {
    int count = 0;
    for (const auto& triangleA : Ta.triangles) {
        for (const auto& triangleB : Tb.triangles) {
            if (areAdjacent(triangleA, triangleB)) {
                count++;
            }
        }
    }
    return count;
}

bool isTriangleAdjacentToCluster(const Triangle& triangle, const Cluster& cluster) {
    for (const auto& t : cluster.triangles) {
        if (areAdjacent(triangle, t)) {
            return true;
        }
    }
    return false;
}

void refineClusters(std::vector<Cluster>& clusters, double tau_N) {
    int iteration = 0;

    SpatialHash spatialHash(0.1);
    preprocessSpatialHash(clusters, spatialHash);

    while (true) {
        std::cout << "Starting iteration " << iteration << " with " << clusters.size() << " clusters." << std::endl;

        auto startTime = std::chrono::high_resolution_clock::now();
        int mergesInThisIteration = 0;
        int splitsInThisIteration = 0;
        size_t pairsProcessed = 0;  // Initialize counter for this iteration

        // Declare a variable to accumulate the total time for the specific loop
        std::chrono::milliseconds totalLoopTime(0);


        // std::map<size_t, std::vector<Triangle>> mergeMap;
        std::map<size_t, std::unordered_set<Triangle>> mergeMap; // Change to unordered_set
        std::set<size_t> processedClusters;
        
        for (size_t i = 0; i < clusters.size(); ++i) {
            
            if (clusters[i].triangles.empty() || processedClusters.count(i)) continue;

            std::unordered_set<int> neighboringClusterIndices = spatialHash.getNeighboringClustersForCluster(clusters[i]);
            std::vector<int> neighboringClusterIndicesVec(neighboringClusterIndices.begin(), neighboringClusterIndices.end());
            auto mergeStartTime = std::chrono::high_resolution_clock::now();
            
            #pragma omp parallel for
            for (size_t idx = 0; idx < neighboringClusterIndicesVec.size(); ++idx) {
                int j = neighboringClusterIndicesVec[idx];
                if (i == j || clusters[j].triangles.empty()) continue;

                bool shouldMerge = false;
                std::vector<Triangle> localMergeList;

                double similarity = ATaTb(clusters[i], clusters[j], spatialHash);
                if (similarity >= tau_N) {
                    for (const auto& triangleA : clusters[i].triangles) {
                        auto potentialNeighborsForA = spatialHash.getPotentialNeighbors(triangleA);
                        for (const auto& triangleB : potentialNeighborsForA) {
                            if (triangleA.isAdjacent(triangleB)) {
                                localMergeList.push_back(triangleA);
                                shouldMerge = true;
                            }
                        }
                    }
                }

                #pragma omp critical
                {
                    if (shouldMerge) {
                        // mergeMap[j].insert(mergeMap[j].end(), localMergeList.begin(), localMergeList.end());
                        mergeMap[j].insert(localMergeList.begin(), localMergeList.end());
                        mergesInThisIteration += localMergeList.size();
                        processedClusters.insert(i);
                        processedClusters.insert(j);
                    }
                }
            }
            auto mergeEndTime = std::chrono::high_resolution_clock::now();  // End timing for the specific loop
            totalLoopTime += std::chrono::duration_cast<std::chrono::milliseconds>(mergeEndTime - mergeStartTime);  // Accumulate the time
        }

        std::cout << "Time spent in the specific loop for iteration " << iteration << ": " << totalLoopTime.count() << " milliseconds." << std::endl;

        std :: cout << "BLAAAAAAAAAAAAAAH" << std::endl;
        
        // Perform merges
        for (auto& [clusterIndex, trianglesToMerge] : mergeMap) {
            clusters[clusterIndex].triangles.insert(clusters[clusterIndex].triangles.end(), trianglesToMerge.begin(), trianglesToMerge.end());
        }

        // Handle Splits
        std::list<Cluster> newClusters;

        #pragma omp parallel
        {
            std::list<Cluster> localNewClusters;

            #pragma omp for nowait
            for (size_t i = 0; i < clusters.size(); ++i) {
                if (mergeMap.count(i)) {
                    Cluster splitCluster;
                    for (const auto& triangle : clusters[i].triangles) {
                        if (mergeMap[i].find(triangle) == mergeMap[i].end()) { // Use find method
                            splitCluster.triangles.push_back(triangle);
                        }
                    }
                    if (!splitCluster.triangles.empty()) {
                        #pragma omp critical
                        {
                            // std::cout << "Triangle from cluster " << i << " not found in merge list. Considering for split." << std::endl;
                            // std::cout << "Splitting cluster " << i << " with " << splitCluster.triangles.size() << " triangles." << std::endl;
                            splitsInThisIteration += splitCluster.triangles.size();
                        }
                        localNewClusters.push_back(splitCluster);
                    }
                }
            }

            #pragma omp critical
            {
                newClusters.splice(newClusters.end(), localNewClusters);
            }
        }

        std::cout << "After " << iteration << " iterations, merged " << mergesInThisIteration << " triangles and split " << splitsInThisIteration << " triangles." << std::endl;
        clusters.insert(clusters.end(), newClusters.begin(), newClusters.end());

        // Remove empty clusters after merges and splits
        clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
            [](const Cluster& cluster) { return cluster.triangles.empty(); }),
        clusters.end());

        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        std::cout << "Iteration " << iteration << " took " << duration << " milliseconds." << std::endl;

        if (mergeMap.empty() && newClusters.empty()) {
            std::cout << "No merges or splits detected. Exiting refinement loop." << std::endl;
            break;
        }

        spatialHash.clear();
        preprocessSpatialHash(clusters, spatialHash);
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
    double tau_N = 0.5;  // Threshold for adjacency similarity
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
    SpatialHash spatialHash(1);  // Assuming 1 is the cell size you want

    // Ensure there are at least 3 triangles for the test
    if (triangles.size() < 3) {
        std::cout << "Not enough triangles for the test." << std::endl;
        return;
    }

    // Use the first 3 triangles for the test
    Triangle t1 = triangles[0];
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
    double tau_S = 0.5;  // Threshold for shape similarity
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
    
    // for (int i = 0; i < 10 && i < allVertices.size(); i++) {
    //     std::cout << "Vertex " << i << ": (" << allVertices[i].x << ", " << allVertices[i].y << ", " << allVertices[i].z << ")" << std::endl;
    // }

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


void saveToPLY(const std::vector<Cluster>& clusters, const std::string& filename) {
    std::ostringstream plyStream;

    // Count total vertices and faces
    size_t totalVertices = clusters.size() * 3;  // 3 vertices per triangle
    size_t totalFaces = clusters.size();  // assuming one triangle per cluster for simplicity

    // PLY header
    plyStream << "ply\n";
    plyStream << "format ascii 1.0\n";
    plyStream << "element vertex " << totalVertices << "\n";
    plyStream << "property float x\n";
    plyStream << "property float y\n";
    plyStream << "property float z\n";
    plyStream << "property uchar red\n";  // Vertex color
    plyStream << "property uchar green\n";  // Vertex color
    plyStream << "property uchar blue\n";  // Vertex color
    plyStream << "element face " << totalFaces << "\n";
    plyStream << "property list uchar int vertex_index\n";
    plyStream << "end_header\n";

    // Write vertices with colors
    size_t vertexIndex = 0;
    for (const auto& cluster : clusters) {
        for (const auto& triangle : cluster.triangles) {
            Vector3D color = triangle.color;
            plyStream << triangle.vec1.x << " " << triangle.vec1.y << " " << triangle.vec1.z << " " 
                      << static_cast<int>(color.x * 255) << " " << static_cast<int>(color.y * 255) << " " << static_cast<int>(color.z * 255) << "\n";
            plyStream << triangle.vec2.x << " " << triangle.vec2.y << " " << triangle.vec2.z << " " 
                      << static_cast<int>(color.x * 255) << " " << static_cast<int>(color.y * 255) << " " << static_cast<int>(color.z * 255) << "\n";
            plyStream << triangle.vec3.x << " " << triangle.vec3.y << " " << triangle.vec3.z << " " 
                      << static_cast<int>(color.x * 255) << " " << static_cast<int>(color.y * 255) << " " << static_cast<int>(color.z * 255) << "\n";
        }
    }

    // Write faces
    for (size_t i = 0; i < totalFaces; ++i) {
        plyStream << "3 " << (i*3) << " " << (i*3 + 1) << " " << (i*3 + 2) << "\n";
    }

    // Save to file
    std::ofstream plyFile(filename + ".ply");
    plyFile << plyStream.str();
    plyFile.close();
}