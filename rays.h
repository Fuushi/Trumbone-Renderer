#ifndef SHADERS_H
#define SHADERS_H

#include <vector>
#include <string>

#include "baseClasses.h"

using namespace std;

struct ShaderInputs {
    int objectID = 0;
    double gloss_diffuse_mix_fac = 0; //0 for neutral (frensel)
    iVec3D material_color = {128,128,128}; //default grey
};

class Mesh {
    public:
    // Empty mesh, no vertices or faces initially
    std::vector<Vec3D> vertices;
    std::vector<iVec3D> faces;

    void update_mesh(const std::vector<Vec3D>& new_vertices, const std::vector<iVec3D>& new_faces) {
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
    iVec3D sky_color = {50,50,60};


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


struct Ray {
    //struct contains complex data related to ray intersections
    Vec3D origin;
    Vec3D euler;
    bool intersect = false; //not implemented?
    double depth = 0.0;
    double frensel = 0.0;
    iVec3D color = {20,25,20}; //20,25,20
    Vec3D intersect_point = {0.0,0.0,0.0};
    Vec3D surface_normal = {0.0,0.0,0.0};
    iVec3D reflection_color;
    Vec3D reflection_vec = {0.0,0.0,0.0};

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
    std::vector<double> vec;
    std::vector<int> color;
    //...
};

struct Lux {
    std::vector<LightingContribution> lighting_contributions;
};

std::vector<int> principled_bdsf(Ray ray, Lux lux, World world, ShaderInputs shader_inputs);

#endif