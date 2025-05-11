#ifndef RAYS_H
#define RAYS_H

#include <vector>
#include <string>
#include <sstream> //for parsing
#include <iostream> //can be removed later
#include <fstream> //for file reading

#include "objects.h"
#include "baseClasses.h"
#include "functions.h" //only for debug

using namespace std;

struct ShaderInputs {
    int objectID = 0;
    int material_index = 0;
    double gloss_diffuse_mix_fac = 0; //0 for neutral (frensel)
    iVec3D material_color = {128,128,128}; //default grey
};

class Mesh {
    public:
    // Empty mesh, no vertices or faces initially
    std::vector<Vec3D> vertices;
    std::vector<iVec3D> faces;

    void update_mesh(const std::vector<Vec3D>& new_vertices, const std::vector<iVec3D>& new_faces) {
        //used to pass in external mesh data
        vertices = new_vertices;
        faces = new_faces;
        return;
    }

    std::string load_file(const std::string& file_path) {
        //loads file as string and returns it
        std::ifstream file(file_path);
        if (!file.is_open()) {
            std::cerr << "Error opening file: " << file_path << std::endl;
            return "";
        }

        // Read the entire file into a string
        std::string file_content((std::istreambuf_iterator<char>(file)),
                                  std::istreambuf_iterator<char>());
        file.close();

        return file_content;
    }

    void parse_obj(const std::string& file_path) {
        // Parse OBJ file and populate vertices and faces
        // This is a placeholder for actual OBJ parsing logic
        // For now, we'll just initialize an empty mesh

        //load file from path as string
        std::string raw_file = load_file(file_path);

        //parse file and extract vertices and faces

        // separate by lines (\n demarkator)
        //...

        //establish stream iterable
        std::istringstream stream(raw_file);
        std::string line;

        //clear mesh before loading
        vertices.clear();
        faces.clear();

        //iterate by line
        while (std::getline(stream, line)) {
            // line dumped into line
            if (line.empty()) continue; //skip empty lines

            std::istringstream  line_stream(line);
            std::string prefix;
            line_stream >> prefix; // read prefix

            // detect known prefixes
            if (prefix == "v") {
                //vertex line, parse and add to vertices

                //extract vertex data
                double x, y, z;
                line_stream >> x >> y >> z;

                //pack to Vec3D
                Vec3D vertex(x, y, z);

                //add to vertex array
                vertices.push_back(vertex);
            }
            else if (prefix == "f") {
                //face line, parse and add to faces

                int v1, v2, v3; //vertex indices cache(s)
                char slash; //slash cache

                // Parse 'f v/v/v' or 'f v//v' or 'f v' to extract vertex indices only
                //lambda function to take f"int/int/int" and return the three ints
                auto parse_vertex = [](const std::string& input) -> std::vector<int> {
                    size_t pos1 = input.find('/');
                    size_t pos2 = input.find('/', pos1 + 1);
            
                    if (pos1 == std::string::npos || pos2 == std::string::npos) {
                        throw std::invalid_argument("Invalid format");
                    }
            
                    int a = std::stoi(input.substr(0, pos1));
                    int b = std::stoi(input.substr(pos1 + 1, pos2 - pos1 - 1));
                    int c = std::stoi(input.substr(pos2 + 1));
            
                    return {a-1, b-1, c-1};
                };
                

                std::string v_str;
                std::vector<int> face_indices = {0,0,0};

                int i=0;
                while (line_stream >> v_str) {
                    
                    //parse
                    std::vector<int> index = parse_vertex(v_str);

                    //put first element into the index of the face indices
                    face_indices[i] = index[0];

                    //other data can be extracted but is discarded for now

                    i++;
                }

                //pack to iVec3D
                iVec3D face_indices_ivec(face_indices);

                //push to faces
                faces.push_back(face_indices_ivec);
            }
        }
        return;
    }

    void load_mesh(const std::string& file_path) {
        // Load mesh from file (e.g., OBJ, STL, etc.)
        // This is a placeholder for actual file loading logic
        // Example: Parse OBJ or other 3D file formats
        // For now, we'll just initialize an empty mesh
        parse_obj(file_path); //interpret file type later
    }

    void inspect() {
        for (int i = 0; i < vertices.size(); i++) {
            std::cout << "Vertex: " << vertices[i].x << " " << vertices[i].y << " " << vertices[i].z << std::endl;
        };
        for (int i = 0; i < faces.size(); i++) {
            std::cout << "Face: " << faces[i].x << " " << faces[i].y << " " << faces[i].z << std::endl;
        };
        //display metadata
        std::cout << "Mesh has: " << vertices.size() << " vertices" << std::endl;
        std::cout << "Mesh has: " << faces.size() << " faces" << std::endl;
    };

    // Constructor that initializes a mesh from a file path
    Mesh(const std::string& file_path) {
        load_mesh(file_path);
    }
    Mesh() = default;
};

//procedural mesh generator, has to be defined after Mesh
class Procedural_meshes {
    public:
    //procedurally generates a planar mesh and assigns it to the Mesh object
    void plane(Mesh& mesh) {
        // Define new vertices and faces
        std::vector<Vec3D> newVertices = {
            {0.0, -6.0, -6.0},
            {0.0, -6.0, 6.0},
            {0.0, 6.0, -6.0},
            {0.0, 6.0, 6.0}
        };


        std::vector<iVec3D> newFaces = {
            {0, 1, 2},//,  // First face
            {2, 1, 3}   // Second face
        };

        mesh.update_mesh(newVertices, newFaces);
        return;
    }

    void neutral_surface(Mesh& mesh) {
        std::vector<Vec3D> vertices = {
            {-1, -1, 0},
            {-1, 1, 0},
            {1, -1, 0},
            {1, 1, 0},

        };

        std::vector<iVec3D> faces = {
            {0, 1, 2},//,  // First face
            {2, 1, 3}   // Second face
        };

        mesh.update_mesh(vertices, faces);
        return;
    }

    void neutral_surface_rotated(Mesh& mesh) {
        std::vector<Vec3D> vertices {
            {-1, 0, -1},
            {-1, 0, 1},
            {1, 0, -1},
            {1, 0, 1}
        };

        std::vector<iVec3D> faces = {
            {0, 1, 2},//,  // First face
            {2, 1, 3}   // Second face
        };
        mesh.update_mesh(vertices, faces);
        return;
    }

    void cube(Mesh& mesh, double size = 1.0) {
        // Half-size to center the cube at the origin
        double halfSize = size / 2.0;

        // Define the vertices of the cube
        std::vector<Vec3D> newVertices = {
            {-halfSize, -halfSize, -halfSize},  // Vertex 0
            {halfSize, -halfSize, -halfSize},   // Vertex 1
            {halfSize, halfSize, -halfSize},    // Vertex 2
            {-halfSize, halfSize, -halfSize},   // Vertex 3
            {-halfSize, -halfSize, halfSize},   // Vertex 4
            {halfSize, -halfSize, halfSize},    // Vertex 5
            {halfSize, halfSize, halfSize},     // Vertex 6
            {-halfSize, halfSize, halfSize}     // Vertex 7
        };

        // Define the faces of the cube using the vertices
        std::vector<iVec3D> newFaces = {
            {0, 1, 2}, {0, 2, 3}, // Front face
            {4, 5, 6}, {4, 6, 7}, // Back face
            {0, 1, 5}, {0, 5, 4}, // Bottom face
            {3, 2, 6}, {3, 6, 7}, // Top face
            {0, 3, 7}, {0, 7, 4}, // Left face
            {1, 2, 6}, {1, 6, 5}  // Right face
        };

        // Update the mesh with the new vertices and faces
        mesh.update_mesh(newVertices, newFaces);
        return;
    }

};

class Element {
    public:
    //Element Container
    Mesh mesh;
    ShaderInputs shaderInputs;


    //object data
    Vec3D pos;

    //update mesh based on pos
    

    //textures
    //...

    //debug method
    void inspect() {
        //for (int i = 0; i < mesh.vertices.size(); i++) {
        //   cout << "Vertex: " << mesh.vertices[i][0] << " " << mesh.vertices[i][1] << " " << mesh.vertices[1][2] << endl;
        //};

        //std::cout << "Mesh has: " << mesh.vertices.size() << " vertices" << endl;
        //std::cout << "Mesh has: " << mesh.faces.size() << " faces" << endl;
    };

    //constructor
    Element(const Mesh& inputMesh) : mesh(inputMesh) {}
};

class Light {
    public:

    std::string type = "sun";
    double brightness = 1.0;
    iVec3D color = {255,255,255};
    Vec3D vec;
};


//define world class
class World {
    public:

    //world contains global render properties
    iVec3D sky_color = {0,50,0}; //default sky color (pre shader)


    //world contains an array of elements
    std::vector<Element> elements;
    std::vector<Light> lights;

    //method to add an element to the world
    void addElement(const Element& element) {
        elements.push_back(element);
    }

    void addLight (const Light& light) {
        lights.push_back(light);
    }

    // Mockup constructor that calls the addElement function
    World() {
        // Create a default mesh and element
        Mesh emptyMesh;
        emptyMesh.faces = {}; // Initialize with empty faces   
        emptyMesh.vertices = {}; // Initialize with empty vertices

        // initialize sky element
        Element skyObj(emptyMesh);
        addElement(skyObj); // Add the default element to the world
    }
};


//define camera class
class Camera {
    public:
    //position
    //std::vector<double> vec3D_c_pos = {3.0, -10.0, 0.0};
    Vec3D vec3D_c_pos = {3.0,-10.0, 0.0};
    //unit vector
    //std::vector<double> vec3D_c_vec = {0.0, 1.0, 0.0};
    Vec3D vec3D_c_vec = {0.0, 1.0, 0.0};

    //camera info
    double fov = 75.0;
    int max_bounce = 1;
    int square_res = 1024;
    iVec2D res = {square_res, square_res};
};

//define state class
class State {
    public:
    World world;
    Camera camera;
};

struct Intersect {
    double depth;
    Vec3D intersect_point;
    Vec3D surface_normal;

    Element& target; //MUST be initialized
    //...
    // Constructor to initialize the reference
    Intersect(Element& element) : target(element) {}
};

struct FastRay {
    bool intersect;
    double depth;
};

struct LightingContribution {
    bool obstructed;
    double depth;
    double brightness;
    Vec3D vec;
    iVec3D color;
    //...
};

struct Lux {
    std::vector<LightingContribution> lighting_contributions;
};

//shader wrapper
std::vector<int> shader_wrapper(Ray ray, Lux lux, World world, ShaderInputs shader_inputs);

std::vector<int> sky_shader(Ray ray);

std::vector<int> geometry_shader(Ray ray, Lux lux, World world, ShaderInputs shader_inputs);

#endif