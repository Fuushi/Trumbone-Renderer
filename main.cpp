#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

//file io declarations
void saveAsPPM(const std::vector<std::vector<std::vector<double>>>& uvs, const std::string& filename) {
    // Get dimensions
    size_t height = uvs.size();
    size_t width = height > 0 ? uvs[0].size() : 0;

    // Open output file
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << "\n";
        return;
    }

    // Write PPM header
    file << "P6\n" << width << " " << height << "\n255\n";

    // Write pixel data
    for (const auto& row : uvs) {
        for (const auto& pixel : row) {
            // Ensure pixel has exactly three components
            //if (pixel.size() != 3) {
            //   std::cerr << "Error: Pixel data must have exactly 3 color channels (R, G, B).\n";
            //   return;
            //}
            // Convert [0.0, 1.0] to [0, 255] and write as binary
            char r = static_cast<char>(pixel[0] * 255);
            char g = static_cast<char>(pixel[1] * 255);
            char b = static_cast<char>(0);
            file.write(&r, 1);
            file.write(&g, 1);
            file.write(&b, 1);
        }
    }

    file.close();
    if (!file) {
        std::cerr << "Error: Failed to properly write to file " << filename << "\n";
    } else {
        std::cout << "Image saved as " << filename << "\n";
    }
}



std::vector<double> matrix_vector_multiplication(std::vector<double> vector1, std::vector<std::vector<double>> vector2) {
    for (int x = 0; x < vector2.size(); x++) {
        vector1[0] = (vector2[x][0] * vector1[0]) + (vector2[x][1] * vector1[0]) + (vector2[x][2] * vector1[0]);
    }
    return vector1;
};

std::vector<double> rotation_matrix_degrees(std::vector<double> input_vector, double fov, std::vector<double> uv) {
    //applies rotation matrix to a vector according to degrees, uv
    double theta_z_deg = (uv[0]*fov) - (fov / 2);
    double theta_x_deg = (uv[1]*fov) - (fov / 2);

    //convert to radians
    double theta_z_rad = theta_z_deg * (3.14159265358979323846264338327950288419716939937510582097494459230781640628 / 180);
    double theta_x_rad = theta_x_deg * (3.14159265358979323846264338327950288419716939937510582097494459230781640628 / 180);

    //define z-axis rotation matrix
    std::vector<std::vector<double>> r_z = {
        {cos(theta_z_rad), 0.0-sin(theta_z_rad), 0.0},
        {sin(theta_z_rad), cos(theta_z_rad), 0.0},
        {0.0, 0.0, 1.0}
    };

    //define x-axis rotation matrix
    std::vector<std::vector<double>> r_x = {
        {1.0, 0.0, 0.0},
        {0.0, cos(theta_x_rad), 0.0-sin(theta_x_rad)},
        {0.0, sin(theta_x_rad), cos(theta_x_rad)}
    };

    //apply matrix(s) via dot product
    //...
    std::vector<double> z_corrected_matrix = matrix_vector_multiplication(input_vector, r_z);
    std::vector<double> x_corrected_matrix = matrix_vector_multiplication(z_corrected_matrix, r_x);


    return x_corrected_matrix;
};





class Mesh {
    public:
    // Empty mesh, no vertices or faces initially
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;

    void update_mesh(std::vector<std::vector<double>> new_vertices, std::vector<std::vector<int>> new_faces) {
        vertices = new_vertices;
        faces = new_faces;
        return;
    }

    // Default constructor that initializes an empty mesh
    Mesh() = default;
};

//procedural mesh generator, has to be defined after Mesh
class Procedural_meshes {
    public:
    //procedurally generates a planar mesh and assigns it to the Mesh object
    void plane(Mesh& mesh) {
        // Define new vertices and faces
        std::vector<std::vector<double>> newVertices = {
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {1.0, 1.0, 0.0},
            {0.0, 1.0, 0.0}
        };

        std::vector<std::vector<int>> newFaces = {
            {0, 1, 2},  // First face
            {0, 2, 3}   // Second face
        };

        mesh.update_mesh(newVertices, newFaces);
        return;
    }
};

class Element {
    public:
    //Element contains
    //mesh
    Mesh mesh;
    

    //textures
    //...

    //debug method
    void inspect() {
        std::cout << "Mesh has: " << mesh.vertices.size() << " vertices" << endl;
        std::cout << "Mesh has: " << mesh.faces.size() << " faces" << endl;
    };

    //constructor
    Element(const Mesh& inputMesh) : mesh(inputMesh) {}
};

//define world class
class World {
    public:
    //world contains an array of elements
    std::vector<Element> elements;

    //method to add an element to the world
    void addElement(const Element& element) {
        elements.push_back(element);
    }

};

//define camera class
class Camera {
    public:
    //position
    std::vector<double> vec3D_c_pos = {-10.0, 0.0, 0.0};
    //unit vector
    std::vector<double> vec3D_c_vec = {10.0, 0.0, 0.0};

    //camera info
    double fov = 90.0;
    int max_bounce = 1;
    std::vector<int> res = {480, 480};

};

//define state class
class State {
    public:
    World world;
    Camera camera;
};

//Rendering
//ray struct for passing into function
struct Ray {
    //DEPRECIATED
};

class OutputBuffer {
    public:
    //Creates an output buffer of empty arrays, pass in as reference
    std::vector<std::vector<std::vector<int>>> rgb_buffer;
    std::vector<std::vector<std::vector<double>>> uvs;

    //sizes buffers accordingly (incorrect buffer sizing fails silently)
    //This can be constructor
    void initializeBuffers(std::vector<int> res) {
        // Ensure the UV buffer is properly initialized
        uvs.resize(res[0]);
        for (auto& row : uvs) {
            row.resize(res[1]);
            for (auto& cell : row) {
                cell.resize(2); // Each cell is a vector with 2 doubles
            }
        }
    };

    void output_buffers_to_file(string fileName) {
        saveAsPPM(uvs, fileName);
    };
};

class Render {
    private:
        State& state;           // Declare member variable
        OutputBuffer& buffer;   // Declare member variable

    public:
    Render(State& state, OutputBuffer& buffer) 
        : state(state), buffer(buffer) { // Initialize members using initializer list
    }


    void generate_uvs() {
        //
        cout << "Generating UVs" << endl;

        int a = 0;
        for (int x = 0; x < state.camera.res[0]; x++) {
            for (int y = 0; y < state.camera.res[1]; y++) {
                //cout << x << " " << y << endl;
                // generate UV and assign to UV buffer
                double uv_x = x / double(state.camera.res[0]);
                double uv_y = y / double(state.camera.res[1]);

                //assign to uv buffer
                buffer.uvs[x][y] = {uv_x, uv_y}; //for some reason, this line stops the next line from printing
                //cout << buffer.uvs[x][y][0] << " " << buffer.uvs[x][y][1] << endl;
            };
        };

    };
    
    void assign_rays() {
        //get camera
        std::vector<double> cam_pos = state.camera.vec3D_c_pos;
        std::vector<double> cam_euler = state.camera.vec3D_c_vec;

        //get UVs
        for (int x = 0; x < state.camera.res[0]; x++) {
            for (int y = 0; y < state.camera.res[1]; y++) {
                //get uvs from buffer
                double uv_x = buffer.uvs[x][y][0];
                double uv_y = buffer.uvs[x][y][1];

                std::vector<double> uv = {uv_x, uv_y};

                //use UV value to calculate the ray
                std::vector<double> ray_euler = rotation_matrix_degrees(cam_euler, state.camera.fov, uv);

                //print ray
                std::cout << ray_euler[0] << " " << ray_euler[1] << " " << ray_euler[2] << endl;
                


            };
        };


    };

    void cast_ray(std::vector<double> ray_origin, std::vector<double> ray_euler) {}; //recursive?


};


//main function
int main() {
    //create mesh generator
    Procedural_meshes generator;
    
    //create state object
    State state;

    //initialize buffer
    OutputBuffer buffer;
    buffer.initializeBuffers(state.camera.res);

    //create renderer
    Render camera(state, buffer);

    //create first element, (plane)
    Mesh myMesh;
    generator.plane(myMesh);
    
    //transfer ownership to Element
    Element element(std::move(myMesh));

    //add Element to the world 
    state.world.addElement(element);

    //inspect element
    state.world.elements[0].inspect();

    
    
    //begin render process

    //generate UVs
    camera.generate_uvs();
    camera.assign_rays();

    //dump to file
    buffer.output_buffers_to_file("./tmp/out.ppm");


    cout << "exit 1" << flush;
    return 1;
}