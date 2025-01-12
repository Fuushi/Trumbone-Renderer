#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <math.h>

using namespace std;

double get_magnitude(std::vector<double> vec) {
    return sqrt((vec[0] * vec[0])+(vec[1] * vec[1])+(vec[2] * vec[2]));
}

std::vector<double> set_magnitude(const std::vector<double>& vec, double target_magnitude) {
    // Compute the magnitude of the input vector
    double magnitude = get_magnitude(vec);

    // Create a new vector and normalize its elements
    std::vector<double> result = vec;
    result[0] = (result[0] / magnitude) * target_magnitude;
    result[1] = (result[1] / magnitude) * target_magnitude;
    result[2] = (result[2] / magnitude) * target_magnitude;

    return result;
}

std::vector<int> convert_vec_double_to_int(std::vector<double> vec) {
    return {static_cast<int>(vec[0]), static_cast<int>(vec[1]), static_cast<int>(vec[2])};
}

//file io declarations
void saveAsPPM(const std::vector<std::vector<std::vector<int>>>& img, const std::string& filename) {
    // Get dimensions
    size_t height = img.size();
    size_t width = height > 0 ? img[0].size() : 0;

    // Open output file
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << "\n";
        return;
    }

    // Write PPM header
    file << "P6\n" << width << " " << height << "\n255\n";

    // Write pixel data
    for (const auto& row : img) {
        for (const auto& pixel : row) {
            // Ensure pixel has exactly three components
            //if (pixel.size() != 3) {
            //   std::cerr << "Error: Pixel data must have exactly 3 color channels (R, G, B).\n";
            //   return;
            //}
            // Convert [0.0, 1.0] to [0, 255] and write as binary
            char r = static_cast<char>(pixel[0]);
            char g = static_cast<char>(pixel[1]);
            char b = static_cast<char>(pixel[2]);
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

void saveAsPPM(const std::vector<std::vector<std::vector<double>>> vec, const std::string filename) {    

    //calculate min max magnitudes
    double min = 0.0;
    double max = 0.0;
    for (int x = 0; x < vec.size(); x++) {
        for (int y = 0; y < vec[0].size(); y++) {
            double magnitude = get_magnitude(vec[x][y]);
            if (magnitude < min) {
                min = magnitude;
            };
            if (magnitude > max) {
                max = magnitude;
            }
            
        }
    };

    //create new array and modulate magnitudes by
    //p=color_bit_range
    //r0 = min value
    //r1 = max value
    //value = (p/(|r0-r1|))*(x-r0) <<< set magnitude to value
    double p = 255.0;
    double r0=min;
    double r1=max;
    double value=0;
    std::vector<std::vector<std::vector<int>>> new_bitmap_array;

    //size bitmap
    new_bitmap_array.resize(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        new_bitmap_array[i].resize(vec[i].size());
        for (size_t j = 0; j < vec[i].size(); ++j) {
            new_bitmap_array[i][j].resize(vec[i][j].size());
        }
    }

    //assign new values
    for (int x = 0; x < vec.size(); x++) {
        for (int y = 0; y < vec[0].size(); y++) {
            //establish target magnitude by function
            value = (p/(r0-r1))*(get_magnitude(vec[x][y])-r0);

            //set magnitude by ref
            new_bitmap_array[x][y] = convert_vec_double_to_int(set_magnitude(vec[x][y], value));
        }
    }

    //call origin
    saveAsPPM(new_bitmap_array, filename);
};

void saveAsPPM_DEPRECIATED(const std::vector<std::vector<std::vector<double>>>& uvs, const std::string& filename) {
    
    //normalize
    
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
            int dim = pixel.size();

            if (dim == 2) {
                char r = static_cast<char>(pixel[0] * 255);
                char g = static_cast<char>(pixel[1] * 255);
                char b = static_cast<char>(0);
                file.write(&r, 1);
                file.write(&g, 1);
                file.write(&b, 1);
            }
            else {
                char r = static_cast<char>(pixel[0]);
                char g = static_cast<char>(pixel[1]);
                char b = static_cast<char>(pixel[2]);
                file.write(&r, 1);
                file.write(&g, 1);
                file.write(&b, 1);
            }
            
        }
    }

    file.close();
    if (!file) {
        std::cerr << "Error: Failed to properly write to file " << filename << "\n";
    } else {
        std::cout << "Image saved as " << filename << "\n";
    }
}

double vector_difference(std::vector<double> vec1, std::vector<double> vec2) {
    //make pythagoras proud
    double sqare_sum = (vec1[0]-vec2[0])*(vec1[0]-vec2[0]) + ((vec1[1]-vec2[1])*(vec1[1]-vec2[1])) + (((vec1[2]-vec2[2])*(vec1[2]-vec2[2]))); 
    return sqrt(sqare_sum);
};

std::vector<double> normalize_vector(const std::vector<double>& vec) {
    // Calculate the magnitude of the vector
    double magnitude = 0.0;
    for (double val : vec) {
        magnitude += val * val;
    }
    magnitude = std::sqrt(magnitude);

    // Check for zero magnitude to avoid division by zero
    if (magnitude == 0.0) {
        throw std::invalid_argument("Cannot normalize a zero-length vector");
    }

    // Create the normalized vector
    std::vector<double> normalized_vec(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        normalized_vec[i] = vec[i] / magnitude;
    }

    return normalized_vec;
}


std::vector<double> vector_cross_product(std::vector<double> vec_a, std::vector<double> vec_b) {
    std::vector<double> cross_product_vector = {
        vec_a[1]*vec_b[2] - vec_a[2]*vec_b[1],
        vec_a[2]*vec_b[0] - vec_a[0]*vec_b[2],
        vec_a[0]*vec_b[1] - vec_a[1]*vec_b[0]
    };
    return cross_product_vector;
};


double vector_dot_product(std::vector<double> vec_a, std::vector<double> vec_b) {
    double product = 0;

    for (int i = 0; i < vec_a.size(); i++) {
        product = product + vec_a[i]*vec_b[i];
    };
    return product;
};


std::vector<double> matrix_subration(std::vector<double> vec_a, std::vector<double> vec_b) {
    std::vector<double> result(vec_a.size());

    for (int i = 0; i < vec_a.size(); i++) {
        result[i] = vec_a[i] - vec_b[i];
    };

    return result;
};

std::vector<double> scalar_multiply(const std::vector<double> vec, double scalar) {
    return {vec[0] * scalar, vec[1] * scalar, vec[2] * scalar};
}

// Function to add two vectors
std::vector<double> vector_add(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    return {vec1[0] + vec2[0], vec1[1] + vec2[1], vec1[2] + vec2[2]};
}

std::vector<double> matrix_vector_multiplication(
    std::vector<double> vector1,
    const std::vector<std::vector<double>>& matrix) 
{
    // Ensure dimensions match for multiplication
    if (matrix.empty() || matrix[0].size() != vector1.size()) {
        throw std::invalid_argument("Matrix and vector dimensions do not match");
    }

    // Resultant vector of size equal to the number of rows in the matrix
    std::vector<double> result(matrix.size(), 0.0);

    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix[i].size(); ++j) {
            result[i] += matrix[i][j] * vector1[j];
        }
    }

    return result;
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

    //normalize
    std::vector<double> normalized_vector = normalize_vector(x_corrected_matrix);
    
    return normalized_vector;
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
            {0.0, -6.0, -6.0},
            {0.0, -6.0, 6.0},
            {0.0, 6.0, -6.0},
            {0.0, 6.0, 6.0}
        };


        std::vector<std::vector<int>> newFaces = {
            {0, 1, 2},//,  // First face
            {2, 1, 3}   // Second face
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
        //for (int i = 0; i < mesh.vertices.size(); i++) {
        //   cout << "Vertex: " << mesh.vertices[i][0] << " " << mesh.vertices[i][1] << " " << mesh.vertices[1][2] << endl;
        //};

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
    std::vector<double> vec3D_c_pos = {3.0, -10.0, 0.0};
    //unit vector
    std::vector<double> vec3D_c_vec = {0.0, 1.0, 0.0};

    //camera info
    double fov = 75.0;
    int max_bounce = 1;
    int square_res = 1024;
    std::vector<int> res = {square_res, square_res};

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
    std::vector<std::vector<std::vector<double>>> rgb_buffer;
    std::vector<std::vector<std::vector<double>>> uvs;
    std::vector<std::vector<std::vector<double>>> ray_euler;

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

        rgb_buffer.resize(res[0]);
        for (auto& row : rgb_buffer) {
            row.resize(res[1]);
            for (auto& cell : row) {
                cell.resize(2); // Each cell is a vector with 2 doubles
            }
        }

        ray_euler.resize(res[0]);
        for (auto& row : ray_euler) {
            row.resize(res[1]);
            for (auto& cell : row) {
                cell.resize(2); // Each cell is a vector with 2 doubles
            }
        }
    };

    void output_buffers_to_file(string fileName) {
        saveAsPPM(uvs, fileName+"uv.ppm");
        saveAsPPM(rgb_buffer, fileName+"mask.ppm");
        saveAsPPM(ray_euler, fileName+"vec.ppm");

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
                //std::cout << ray_euler[0] << " " << ray_euler[1] << " " << ray_euler[2] << endl;
                
                //cast ray via function call

                //create rgb pixel
                std::vector<double> pixel = cast_ray(cam_pos, ray_euler);

                //assign to buffer
                buffer.rgb_buffer[x][y] = pixel;
                buffer.ray_euler[x][y] = ray_euler;


            };
        };


    };

    std::vector<double> cast_ray(std::vector<double> ray_origin, std::vector<double> ray_euler) {

        //define intersects array
        std::vector<std::vector<double>> intersects;
        

        //iterate through world assets
        for (int i = 0; i < state.world.elements.size(); i++) {
            
            //get asset elements
            Element element = state.world.elements[i];

            //iterate through faces
            std::vector<int> face;
            std::vector<std::vector<double>> vertices;
            for (int j = 0; j < element.mesh.faces.size(); j++) {

                //extract face
                face = element.mesh.faces[j];

                //extract vertices
                vertices = {
                    element.mesh.vertices[face[0]],
                    element.mesh.vertices[face[1]],
                    element.mesh.vertices[face[2]]
                };

                //define edges
                std::vector<double> edge1 = matrix_subration(vertices[1], vertices[0]); //v2-v1
                std::vector<double> edge2 = matrix_subration(vertices[2], vertices[0]); //v3-v1

                //compute normal of the 2 edges
                //... (cross(e1,e2))
                std::vector<double> surface_normal = vector_cross_product(edge1, edge2);

                //compute determinant of edge 1 and 2
                //...

                //pvec = cross(ray_vec, e2)
                std::vector<double> pvec = vector_cross_product(ray_euler, edge2);

                //det = dot(e1, pvec)
                double determinant = vector_dot_product(edge1, pvec);

                //ray parrallel optimization
                if (determinant == 0) {
                    //cout << "Ray Parralell";
                    continue;
                };

                //compute inverse determinant
                double inverse_determinant =  1 / determinant;

                //compute u in normal space
                std::vector<double> tvec = matrix_subration(ray_origin, vertices[0]);
                double u = vector_dot_product(tvec, pvec) * inverse_determinant;
                if ((u < 0) || (u > 1)) {
                    //no intersection
                    continue;
                }

                //compute v in normal space
                std::vector<double> qvec = vector_cross_product(tvec, edge1);
                double v = vector_dot_product(ray_euler, qvec) * inverse_determinant;
                if ((v < 0) || ((u + v) > 1)) {
                    //no intersection
                    continue;
                }

                //compute t in normal space
                double t = vector_dot_product(edge2, qvec) * inverse_determinant;
                if (t < 0) {
                    //no intersect
                    continue;
                }

                //intersect

                //compute depth
                //double depth = vector_difference(ray_origin, {u,v,t});
                //cout << "Intersect";
                

                //... compute intersect data
                std::vector<double> mp = scalar_multiply(ray_euler, t);
                std::vector<double> intersection_point = vector_add(mp, ray_origin);

                return intersection_point;

            };
        }
        return {0,0,0};   
    }; //recursive?


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
    buffer.output_buffers_to_file("./tmp/");


    return 1;

    for (double i = 0; i < 60*10; i=i+1.0) {
                //generate UVs
        camera.generate_uvs();
        std::vector<double> rotation_uv_lmfao = {0.5, 0.6};
        state.camera.vec3D_c_vec = rotation_matrix_degrees(
            state.camera.vec3D_c_vec, 75.0, rotation_uv_lmfao
        );
        camera.assign_rays();

        //dump to file
        buffer.output_buffers_to_file("./tmp/");
    }

    cout << "exit 1" << flush;
    return 1;
}