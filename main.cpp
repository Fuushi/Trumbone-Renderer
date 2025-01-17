#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <math.h>

//functions.h
#include "functions.h"
#include "rays.h"

using namespace std;


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


class OutputBuffer {
    public:
    //Creates an output buffer of empty arrays, pass in as reference
    std::vector<std::vector<std::vector<int>>> rgb_buffer;
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
        saveAsPPM(rgb_buffer, fileName+"int.ppm");
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

                //create rgb pixel
                Ray ray = cast_ray(cam_pos, ray_euler);


                //assign to buffer
                buffer.rgb_buffer[x][y] = ray.color;
                buffer.ray_euler[x][y] = ray.reflection_vec; //composite pass


            };
        };


    };

    Ray cast_ray(std::vector<double> ray_origin, std::vector<double> ray_euler, int recursion_depth=0) {

        //decide recursion
        if (recursion_depth > state.camera.max_bounce+1) {
            Ray ray;
            return ray;
        }

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

                //... compute intersect data
                std::vector<double> mp = scalar_multiply(ray_euler, t);
                std::vector<double> intersection_point = vector_add(mp, ray_origin);

                //calculate depth
                double depth = vector_difference(intersection_point, ray_origin);

                //calculate frensel
                double frensel_deg = 180.0 - get_vectors_angle(surface_normal, ray_euler);

                //calculate normal
                //precomputed at surface_normal

                //calculate reflection
                std::vector<double> reflection_vector = reflect_vector_normal(ray_euler, surface_normal);

                //cast next bounce
                double epsilon = 1e-9; //advance ray forward by some epsilon
                std::vector<double> advanced_ray = vector_add(ray_origin, set_magnitude(reflection_vector, epsilon));

                Ray reflect_ray = cast_ray(
                    advanced_ray, 
                    reflection_vector, 
                    recursion_depth+1
                );

                //pack struct
                Ray ray;
                ray.origin = ray_origin;
                ray.euler = ray_euler;
                ray.intersect=true;
                ray.intersect_point = intersection_point;
                ray.surface_normal = surface_normal;
                ray.depth = depth;
                ray.frensel = frensel_deg;
                ray.reflection_vec = reflection_vector;
                ray.reflection_color = reflect_ray.color;

                //shader
                std::vector<int> pixel = principled_bdsf(ray, state.world);

                ray.color = pixel;

                return ray;

            };
        }
        Ray ray; //empty ray
        return ray;
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

    //create lights
    Light light;
    light.vec={-1.0,-1.0,-1.0};

    //add light to scene
    state.world.addLight(light);

    
    
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