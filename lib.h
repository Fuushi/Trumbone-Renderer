#ifndef LIB_H
#define LIB_H

#include <vector>
#include <string>
#include <iostream>

//functions.h
#include "functions.h"
#include "rays.h"

using namespace std;


void saveAsPPM(const std::vector<std::vector<std::vector<int>>>& img, const std::string& filename);
void saveAsPPM(const std::vector<std::vector<std::vector<double>>> vec, const std::string filename);


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
        std::vector<Intersect> intersects;
        

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
                std::vector<double> intersection_point = vector_add(
                    scalar_multiply(ray_euler, t), ray_origin
                );

                //calculate depth
                double depth = vector_difference(intersection_point, ray_origin);

                //add intersect to array
                Intersect intersect;
                intersect.depth = depth;
                intersect.intersect_point = intersection_point;
                intersect.surface_normal = surface_normal; 

                intersects.push_back(intersect);

            };
        }



        if (intersects.size() == 0) {
            Ray ray;
            return ray;
        }
        //else

        //select closest intersect
        // Select the closest intersect
        Intersect closestIntersect;
        closestIntersect.depth = std::numeric_limits<double>::max();  // Initialize with max double value

        for (int i = 0; i < intersects.size(); i++) {
            if (intersects[i].depth < closestIntersect.depth) {  // Check for the closest (smallest depth)
                closestIntersect = intersects[i];
            }
        }


        //calculate frensel
        double frensel_deg = 180.0 - get_vectors_angle(closestIntersect.surface_normal, ray_euler);

        //calculate normal
        //precomputed at surface_normal

        //compute illumination...

        //calculate reflection
        std::vector<double> reflection_vector = reflect_vector_normal(ray_euler, closestIntersect.surface_normal);

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
        ray.intersect_point = closestIntersect.intersect_point;
        ray.surface_normal = closestIntersect.surface_normal;
        ray.depth = closestIntersect.depth;
        ray.frensel = frensel_deg;
        ray.reflection_vec = reflection_vector;
        ray.reflection_color = reflect_ray.color;

        //shader
        std::vector<int> pixel = principled_bdsf(ray, state.world);

        ray.color = pixel;

        return ray;
    }; //recursive?


};

#endif