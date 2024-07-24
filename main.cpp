#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <limits>
#include <cmath>

struct Vec2i{
    int x, y;
    Vec2i(int x, int y) : x(x), y(y) {}
    Vec2i() : x(0), y(0) {}
};

struct Vec2f{
    float x, y;
    Vec2f(float x, float y) : x(x), y(y) {}
    Vec2f() : x(0), y(0) {}
};

struct Vec3f{
    float x, y, z;
    Vec3f(float x, float y, float z) : x(x), y(y), z(z) {}
    Vec3f() : x(0), y(0), z(0) {}
    Vec3f operator^(const Vec3f& v) const{
        return Vec3f(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    Vec3f normalize(){
        float norm = sqrt(x * x + y * y + z * z);
        x /= norm;
        y /= norm;
        z /= norm;
        return *this;
    }
    float operator*(const Vec3f& v) const{
        return x * v.x + y * v.y + z * v.z;
    }
    Vec3f operator-(const Vec3f& v) const{
        return Vec3f(x - v.x, y - v.y, z - v.z);
    }
    float& operator[](int i) {
    if (i == 0) return x;
    if (i == 1) return y;
    if (i == 2) return z;
    throw std::out_of_range("Index out of range for Vec3f");
}

};

struct Matrix4x4 {
    float m[4][4] = {0};

    static Matrix4x4 identity() {
        Matrix4x4 result;
        for (int i = 0; i < 4; ++i) {
            result.m[i][i] = 1.0f;
        }
        return result;
    }

    Matrix4x4 operator*(const Matrix4x4& other) const {
        Matrix4x4 result;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    result.m[i][j] += m[i][k] * other.m[k][j];
                }
            }
        }
        return result;
    }

    Vec3f operator*(const Vec3f& v) const {
        Vec3f result;
        float w;
        result.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3];
        result.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3];
        result.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3];
        w = m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3];

        if (w != 1 && w != 0) {
            result.x /= w;
            result.y /= w;
            result.z /= w;
        }

        return result;
    }
};




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
    void normalize() {
        // 初始化每个轴的最小和最大值
        float minX = std::numeric_limits<float>::max();
        float maxX = std::numeric_limits<float>::lowest();
        float minY = std::numeric_limits<float>::max();
        float maxY = std::numeric_limits<float>::lowest();
        float minZ = std::numeric_limits<float>::max();
        float maxZ = std::numeric_limits<float>::lowest();

        // 找到每个轴的最小和最大值
        for (const auto& facet : facets) {
            for (int i = 0; i < 3; ++i) {
                minX = std::min(minX, facet.vertices[i].x);
                maxX = std::max(maxX, facet.vertices[i].x);
                minY = std::min(minY, facet.vertices[i].y);
                maxY = std::max(maxY, facet.vertices[i].y);
                minZ = std::min(minZ, facet.vertices[i].z);
                maxZ = std::max(maxZ, facet.vertices[i].z);
            }
        }
        // print min max xyz
        // std::cout << "minX: " << minX << " maxX: " << maxX << std::endl;
        // std::cout << "minY: " << minY << " maxY: " << maxY << std::endl;
        // std::cout << "minZ: " << minZ << " maxZ: " << maxZ << std::endl;

        // 计算每个轴的范围和偏移量
        float rangeX = maxX - minX;
        float scaleX = 2 / rangeX;
        float biasX = (maxX + minX) / 2;
        float rangeY = maxY - minY;
        float scaleY = 2 / rangeY;
        float biasY = (maxY + minY) / 2;
        float rangeZ = maxZ - minZ;
        float scaleZ = 2 / rangeZ;
        float biasZ = (maxZ + minZ) / 2;

        //print scale
        // std::cout << "scaleX: " << scaleX << " biasX: " << biasX << std::endl;
        // std::cout << "scaleY: " << scaleY << " biasY: " << biasY << std::endl;
        // std::cout << "scaleZ: " << scaleZ << " biasZ: " << biasZ << std::endl;

        float scale = std::min(scaleX, std::min(scaleY, scaleZ));
        // 应用归一化
        for (auto& facet : facets) {
            for (int i = 0; i < 3; ++i) {
                facet.vertices[i].x = (facet.vertices[i].x - biasX) * scale;
                facet.vertices[i].y = (facet.vertices[i].y - biasY) * scale;
                facet.vertices[i].z = (facet.vertices[i].z - biasZ) * scale;
            }
        }
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

struct Color{
    uint8_t r, g, b;
    Color(uint8_t r, uint8_t g, uint8_t b) : r(r), g(g), b(b) {}
    Color() : r(0), g(0), b(0) {}
};

class FrameBuffer{
private:
    std::vector<Color> buffer;

public:
    int width, height;
    
    FrameBuffer(int width, int height) : width(width), height(height){
        buffer.resize(width * height);
    }

    void setPixel(int x, int y, const Color& color){
        if(x >= 0 && x < width && y >= 0 && y < height){
            buffer[y * width + x] = color;
        }
    }

    void writePPM(const std::string& filename){
        std::ofstream out(filename, std::ios::binary);
        out << "P6\n" << width << " " << height << "\n255\n";
        for(const auto& color : buffer){
            out.write(reinterpret_cast<const char*>(&color), sizeof(Color));
        }
        out.close();
    }
};

// Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P){
//     Vec3f u = Vec3f(C.x - A.x, B.x - A.x, A.x - P.x) ^ Vec3f(C.y - A.y, B.y - A.y, A.y - P.y);
//     if(std::abs(u.z) > 1e-2) return Vec3f(-1, 1, 1);
//     return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z); 
// }

Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
    // Create vectors from the triangle vertices to the point P
    Vec3f AB = B - A;
    Vec3f AC = C - A;
    Vec3f PA = A - P;

    // Calculate the normal vector using the cross product of vectors AB and AC
    Vec3f n = AB ^ AC;

    // Compute cross products of vectors to P and triangle edges
    Vec3f u = Vec3f(AC.y * PA.z - AC.z * PA.y, AC.z * PA.x - AC.x * PA.z, AC.x * PA.y - AC.y * PA.x);
    Vec3f v = Vec3f(AB.y * PA.z - AB.z * PA.y, AB.z * PA.x - AB.x * PA.z, AB.x * PA.y - AB.y * PA.x);

    // The area of the parallelogram spanned by AB and AC
    float denom = n * n; // Dot product with itself gives the area squared
    if (std::abs(denom) < 1e-2) {
        // Degenerate triangle check (area close to zero)
        //std::cout << "Degenerate triangle" << std::endl;
        return Vec3f(-1, 1, 1); // Use a special value to indicate a degenerate case
    }

    // The area ratios, subtracted from 1 due to how u and v are computed relative to point P
    float uArea = (n * u) / denom; // Project u onto n and divide by total area
    float vArea = (n * v) / denom; // Project v onto n and divide by total area
    float wArea = 1.0f - uArea - vArea;
    //std::cout << "uArea: " << uArea << " vArea: " << vArea << " wArea: " << wArea << std::endl;
    return Vec3f(wArea, uArea, vArea);
}

class ZBuffer{
private:
    float* buffer;
    int width, height;
public:
    ZBuffer(int width, int height) : width(width), height(height){
        buffer = new float[width * height];
        for(int i = 0; i < width * height; i++){
            buffer[i] = std::numeric_limits<float>::lowest();
        }
    }

    ~ZBuffer(){
        delete[] buffer;
    }

    bool set(int x, int y, float z){
        if(x >= 0 && x < width && y >= 0 && y < height){
            if(z > buffer[y * width + x]){
                buffer[y * width + x] = z;
                return true;
            }  
        }
        return false;
    }

    float get(int x, int y){
        if(x >= 0 && x < width && y >= 0 && y < height){
            return buffer[y * width + x];
        }
        return std::numeric_limits<float>::lowest();
    }
};

struct Camera{
    Vec3f position;
    Vec3f direction;
    Vec3f up;
    float angle;

    Camera(Vec3f position, Vec3f direction, Vec3f up, float angle) : position(position), direction(direction), up(up), angle(angle) {}

    Matrix4x4 getViewMatrix() {
        Vec3f zAxis = direction.normalize();
        Vec3f xAxis = (up ^ zAxis).normalize();
        Vec3f yAxis = zAxis ^ xAxis;

        Matrix4x4 view = Matrix4x4::identity();
        view.m[0][0] = xAxis.x; view.m[0][1] = xAxis.y; view.m[0][2] = xAxis.z; view.m[0][3] = -(xAxis * position);
        view.m[1][0] = yAxis.x; view.m[1][1] = yAxis.y; view.m[1][2] = yAxis.z; view.m[1][3] = -(yAxis * position);
        view.m[2][0] = zAxis.x; view.m[2][1] = zAxis.y; view.m[2][2] = zAxis.z; view.m[2][3] = -(zAxis * position);
        return view;
    }

    Matrix4x4 getProjectionMatrix(float aspectRatio) const {
        Matrix4x4 projection = Matrix4x4::identity();
        float near = 0.1f;
        float far = 1000.0f;
        float fov = 1.0f / tan(angle * 0.5 * M_PI / 180);
        projection.m[0][0] = fov / aspectRatio;
        projection.m[1][1] = fov;
        projection.m[2][2] = -(far + near) / (far - near);
        projection.m[2][3] = -2 * far * near / (far - near);
        projection.m[3][2] = -1;
        projection.m[3][3] = 0;

        return projection;
    }


};

int main() {
    STLReader reader("CircleofWillis_clipped.stl");
    reader.normalize();
    const auto& facets = reader.getFacets();
    Camera camera(Vec3f(0, 0, 1), Vec3f(0, 0, -1), Vec3f(0, 1, 0), 90);

    int width = 800;
    int height = 800;

    Matrix4x4 view = camera.getViewMatrix();
    Matrix4x4 projection = camera.getProjectionMatrix((float)width / height);
    Matrix4x4 transform = projection * view;

    FrameBuffer framebuffer(width, height);
    ZBuffer zbuffer(width, height);
    // Draw the facets
    for (const auto& facet : facets) {
        Vec3f screen_coords[3];
        Vec3f camera_coords[3];
        Vec3f world_coords[3];
        for (int i = 0; i < 3; ++i) {
            // 模型空间 -> 世界空间
            world_coords[i] = Vec3f(facet.vertices[i].x, facet.vertices[i].y, facet.vertices[i].z);
            // print world_coords
            //std::cout << "world_coords[" << i << "]: " << world_coords[i].x << " " << world_coords[i].y << " " << world_coords[i].z << std::endl;
            camera_coords[i] = view * world_coords[i];

            screen_coords[i] = transform * world_coords[i];
            screen_coords[i].x = (screen_coords[i].x + 1) * width / 2;
            screen_coords[i].y = (screen_coords[i].y + 1) * height / 2;
            // print screen_coords
            //std::cout << "screen_coords[" << i << "]: " << screen_coords[i].x << " " << screen_coords[i].y << " " << screen_coords[i].z << std::endl;
        }
        Vec3f n = (camera_coords[2] - camera_coords[0]) ^ (camera_coords[1] - camera_coords[0]);
        n.normalize();
        float intensity = n * Vec3f(0, 0, -1);
        // Draw the edges
        Color color(128*intensity, 0, 0);
        Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        Vec2f bboxmax(std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest());
        for(int i = 0; i < 3; i++){
            bboxmin.x = std::min(bboxmin.x, screen_coords[i].x);
            bboxmin.y = std::min(bboxmin.y, screen_coords[i].y);
            bboxmax.x = std::max(bboxmax.x, screen_coords[i].x);
            bboxmax.y = std::max(bboxmax.y, screen_coords[i].y);
        }
        Vec3f P;
        for(P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
            for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
                Vec3f bc_screen = barycentric(screen_coords[0], screen_coords[1], screen_coords[2], P);
                if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
                P.z = 0;
                for(int i = 0; i < 3; i++){
                    P.z += camera_coords[i].z * bc_screen[i];
                }
                int x = P.x;
                int y = P.y;
                if(zbuffer.set(x, y, P.z)){
                    // std::cout << "x: " << x << std::endl;
                    // std::cout << "P.z: " << P.z << std::endl;
                    framebuffer.setPixel(x, y, color);
                }
            }
        }
    }
    framebuffer.writePPM("output.ppm");

    return 0;
}
