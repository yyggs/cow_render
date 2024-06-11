#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

struct Vertex {
    float x, y, z;
};

struct Facet {
    Vertex normal;
    Vertex vertices[3];
};

class STLReader {
public:
    STLReader(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        std::string line;
        std::getline(file, line); // Read the first line to skip the solid name

        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string token;
            ss >> token;

            if (token == "facet") {
                Facet facet;
                ss >> token; // skip "normal"
                ss >> facet.normal.x >> facet.normal.y >> facet.normal.z;
                readFacetVertices(file, facet);
                facets.push_back(facet);
            } else if (token == "endsolid") {
                break;
            }
        }

        file.close();
    }

    const std::vector<Facet>& getFacets() const {
        return facets;
    }

private:
    std::vector<Facet> facets;

    void readFacetVertices(std::ifstream& file, Facet& facet) {
        std::string line;
        std::getline(file, line); // Skip "outer loop"
        for (int i = 0; i < 3; ++i) {
            std::getline(file, line);
            std::stringstream ss(line);
            std::string token;
            ss >> token; // skip "vertex"
            ss >> facet.vertices[i].x >> facet.vertices[i].y >> facet.vertices[i].z;
        }
        std::getline(file, line); // Skip "endloop"
        std::getline(file, line); // Skip "endfacet"
    }
};

int main() {
    STLReader reader("CircleofWillis_clipped.stl");
    const auto& facets = reader.getFacets();

    for (const auto& facet : facets) {
        std::cout << "Normal: " << facet.normal.x << ", " << facet.normal.y << ", " << facet.normal.z << std::endl;
        for (int i = 0; i < 3; ++i) {
            std::cout << "Vertex " << i << ": " 
                      << facet.vertices[i].x << ", " 
                      << facet.vertices[i].y << ", " 
                      << facet.vertices[i].z << std::endl;
        }
    }

    return 0;
}
