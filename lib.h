#ifndef LIB_H
#define LIB_H

#include <vector>
#include <string>
#include <iostream>
#include <limits>

//functions.h
#include "functions.h"
#include "rays.h"
#include "mathHelper.h"
#include "baseClasses.h"

using namespace std;

//function to dump bitmap img to disc
void saveAsPPM(const std::vector<std::vector<std::vector<int>>>& img, const std::string& filename);
//vector bitmap overflow
void saveAsPPM(const std::vector<std::vector<std::vector<double>>> vec, const std::string filename);

//Object to store bitmap img buffers
class OutputBuffer {
    public:
    //Creates an output buffer of empty arrays, pass in as reference
    std::vector<std::vector<std::vector<int>>> rgb_buffer;
    std::vector<std::vector<std::vector<int>>> int_array_2;
    std::vector<std::vector<std::vector<double>>> uvs;
    std::vector<std::vector<std::vector<double>>> ray_euler;

    //sizes buffers accordingly (incorrect buffer sizing fails silently)
    //This can be constructor
    void initializeBuffers(iVec2D res) {
        // Ensure the UV buffer is properly initialized
        uvs.resize(res.x);
        for (auto& row : uvs) {
            row.resize(res.y);
            for (auto& cell : row) {
                cell.resize(2); // Each cell is a vector with 2 doubles
            }
        }

        rgb_buffer.resize(res.x);
        for (auto& row : rgb_buffer) {
            row.resize(res.y);
            for (auto& cell : row) {
                cell.resize(2); // Each cell is a vector with 2 doubles
            }
        }
        
        int_array_2.resize(res.x);
        for (auto& row : int_array_2) {
            row.resize(res.y);
            for (auto& cell : row) {
                cell.resize(2); // Each cell is a vector with 2 doubles
            }
        }

        ray_euler.resize(res.x);
        for (auto& row : ray_euler) {
            row.resize(res.y);
            for (auto& cell : row) {
                cell.resize(2); // Each cell is a vector with 2 doubles
            }
        }
    };

    void output_buffers_to_file(string fileName, bool do_meta_passes=true) {
        saveAsPPM(rgb_buffer, fileName+"composite.ppm");
        if (!do_meta_passes) {return;};

        saveAsPPM(uvs, fileName+"uv.ppm");
        saveAsPPM(ray_euler, fileName+"vec.ppm");
        saveAsPPM(int_array_2, fileName+"int_array_2.ppm");

    };
};

class Animator;

//Renderer Object
class Render {
    private:
        State& state;           // Declare member variable
        OutputBuffer& buffer;   // Declare member variable

    public:
    Render(State& state, OutputBuffer& buffer) 
        : state(state), buffer(buffer) { // Initialize members using initializer list
    }
    
    friend class Animator;  // Grant full access to Animator

    //uv render pass
    void generate_uvs() {
        //
        cout << "Generating UVs..." << endl;

        int a = 0;
        for (int x = 0; x < state.camera.res.x; x++) {
            for (int y = 0; y < state.camera.res.y; y++) {
                //cout << x << " " << y << endl;
                // generate UV and assign to UV buffer
                double uv_x = x / double(state.camera.res.x);
                double uv_y = y / double(state.camera.res.y);

                //assign to uv buffer
                buffer.uvs[x][y] = {uv_x, uv_y}; //for some reason, this line stops the next line from printing
                //cout << buffer.uvs[x][y][0] << " " << buffer.uvs[x][y][1] << endl;
            };
        };

    };
    
    //raycasting render pass
    void assign_rays() {
        std::cout << "Casting Rays..." << std::endl;
        //get camera
        //std::vector<double> cam_pos = state.camera.vec3D_c_pos;
        //std::vector<double> cam_euler = state.camera.vec3D_c_vec;
        Vec3D cam_pos = state.camera.vec3D_c_pos;
        Vec3D cam_euler = state.camera.vec3D_c_vec;


        //get UVs
        for (int x = 0; x < state.camera.res.x; x++) {
            for (int y = 0; y < state.camera.res.y; y++) {
                //get uvs from buffer
                double uv_x = buffer.uvs[x][y][0];
                double uv_y = buffer.uvs[x][y][1];

                std::vector<double> uv = {uv_x, uv_y};

                //use UV value to calculate the ray
                Vec3D ray_euler = calculate_ray_heading(cam_euler, state.camera.fov, uv);

                //create rgb pixel
                Ray ray = cast_ray(cam_pos, ray_euler);

                //assign to buffer
                buffer.rgb_buffer[x][y] = ray.color.to_vec();
                buffer.ray_euler[x][y] = ray.euler.to_double(); //reflection vector

                //
                if (!ray.reflection_color.empty) {
                    buffer.int_array_2[x][y] = ray.reflection_color.to_vec(); //reflection color
                } else {
                    buffer.int_array_2[x][y] = ray.color.to_vec(); //sky color
                }

            };
        };


    };

    //wrapper to combine render passes into a single call
    void render() {
        generate_uvs();
        assign_rays();
    };

    FastRay fast_ray_cast(const Vec3D& origin, const Vec3D& vec) {
        //so long as i never have to touch this code again, it is maintainable
        //returns first detected intersect
        //iterate through meshes


        for (int i = 0; i < state.world.elements.size(); i++) {
            Element element = state.world.elements[i];
            std::vector<Vec3D> vertices;
            iVec3D face;
            //iterate through faces
            for (int j = 0; j < element.mesh.faces.size(); j++) {

                face = element.mesh.faces[j]; //extract mesh data
                vertices = {element.mesh.vertices[face.x],element.mesh.vertices[face.y],element.mesh.vertices[face.z]};

                //define edges
                Vec3D edge1 = vertices[1] - vertices[0]; //v2-v1
                Vec3D edge2 = vertices[2] - vertices[0]; //v3-v1

                //math... elaborated in slow function
                Vec3D pvec = vec.cross(edge2); //pvec = cross(ray_vec, e2)
                double determinant = edge1.dot(pvec); //det = dot(e1, pvec)
                if (determinant == 0) {continue;}; //ray parrallel optimization
                double inverse_determinant =  1 / determinant;
                Vec3D tvec = Vec3D(origin)-vertices[0]; //tvec = ray_origin - v1
                double u = tvec.dot(pvec) * inverse_determinant; //compute u in normal space
                if ((u < 0) || (u > 1)) {continue;}
                Vec3D qvec = tvec.cross(edge1); //qvec = cross(tvec, e1)
                double v = vec.dot(qvec) * inverse_determinant; //compute v in normal space
                if ((v < 0) || ((u + v) > 1)) {continue;}
                double t = edge2.dot(qvec) * inverse_determinant; //compute t in normal space
                if (t < 0) {continue;}
                //compute intersection
                Vec3D intersection_point = origin + (vec % t);//vector_add(scalar_multiply(vec.to_double(), t), origin);

                //calculate depth
                double depth = vector_difference(intersection_point, origin);
                FastRay ray;
                ray.intersect=true;
                ray.depth=depth;
                return ray;
            };
        };
        //no intersects
        FastRay ray;
        ray.intersect=false;
        ray.depth=-1.0;
        return ray;
    };

    LightingContribution cast_light_ray(Light light, const Vec3D& ray_origin, const Vec3D& ray_vec) {
        //calculate ray
        FastRay light_ray = fast_ray_cast(ray_origin, ray_vec);

        //initialize return object
        LightingContribution contribution;

        //pass in ray attributes
        contribution.obstructed = light_ray.intersect;
        contribution.depth = light_ray.depth;

        //pass in light object attributes
        contribution.brightness = light.brightness;
        contribution.vec = light.vec;
        contribution.color = light.color;

        //return LightingContribution object
        return contribution;
    };

    Lux calculate_luminance(const Vec3D& origin) {
        //...
        std::vector<LightingContribution> lighting_contributions;
        for (int i = 0; i < state.world.lights.size(); i++) {
            //

            //light vector is inverted with scalar multiply
            //std::vector<double> inverse_vector = scalar_multiply(state.world.lights[i].vec.to_double(), -1.0);
            Vec3D inverse_vec = state.world.lights[i].vec%-1.0;

            //advance ray by epsilon
            double epsilon = 1e-9; //advance ray forward by some epsilon
            Vec3D advanced_ray = origin + (inverse_vec++ % epsilon);
            //std::vector<double> advanced_ray = vector_add(origin, set_magnitude(inverse_vector, epsilon));

            //compute lighting
            lighting_contributions.push_back(
                cast_light_ray(
                    state.world.lights[i], advanced_ray, inverse_vec
            ));
        };

        // create stuct and assign array by reference
        Lux lux;
        lux.lighting_contributions=lighting_contributions;

        //return struct
        return lux;
    };



    //ray cast function call, returns Ray struct
    Ray cast_ray(const Vec3D& ray_origin, const Vec3D& ray_euler, int recursion_depth=0) {

        //decide recursion
        if (recursion_depth > state.camera.max_bounce+1) {
            Ray ray;
            ray.euler = ray_euler; //set ray euler to the current ray euler (needed for sky shader)
            return ray;
        }

        //define intersects array
        std::vector<Intersect> intersects;
        

        //iterate through world assets
        for (int i = 0; i < state.world.elements.size(); i++) {
            
            //get asset elements
            Element& element = state.world.elements[i];

            //iterate through faces
            iVec3D face;
            std::vector<Vec3D> vertices;
            for (int j = 0; j < element.mesh.faces.size(); j++) {

                //extract face
                face = element.mesh.faces[j];

                //extract vertices
                vertices = {
                    element.pos+element.mesh.vertices[face.x],
                    element.pos+element.mesh.vertices[face.y],
                    element.pos+element.mesh.vertices[face.z]
                };

                //define edges
                Vec3D edge1(vertices[1]-vertices[0]); //v2-v1
                Vec3D edge2(vertices[2]-vertices[0]); //v3-v1

                //compute normal of the 2 edges (cross(e1,e2))
                Vec3D surface_normal(edge1.cross(edge2));

                //pvec = cross(ray_vec, e2)
                Vec3D pvec(ray_euler.cross(edge2));

                //det = dot(e1, pvec)
                double determinant = edge1.dot(pvec);

                //ray parrallel optimization
                if (determinant == 0) {
                    //cout << "Ray Parralell";
                    continue;
                };

                //compute inverse determinant
                double inverse_determinant =  1 / determinant;

                //compute u in normal space
                Vec3D tvec(ray_origin-vertices[0]);
                double u = tvec.dot(pvec) * inverse_determinant;
                if ((u < 0) || (u > 1)) {
                    //no intersection
                    continue;
                }

                //compute v in normal space
                Vec3D qvec(tvec.cross(edge1));
                double v = ray_euler.dot(qvec) * inverse_determinant;
                if ((v < 0) || ((u + v) > 1)) {
                    //no intersection
                    continue;
                }

                //compute t in normal space
                double t = edge2.dot(qvec) * inverse_determinant;
                if (t < 0) {
                    //no intersect
                    continue;
                }

                //... compute intersect data
                Vec3D intersection_point((ray_euler%t)+ray_origin);

                //calculate depth
                double depth = vector_difference(intersection_point, ray_origin);

                //add intersect to array
                Intersect intersect(element);
                intersect.depth = depth;
                intersect.intersect_point = intersection_point;
                intersect.surface_normal = surface_normal; 

                intersects.push_back(intersect);

            };
        }


        // No intersect optimization (skips shader and goes straight to sky to save performance)
        if (intersects.size() == 0) {
            Ray ray;

            //add ray attributes for shaders
            ray.euler = ray_euler; //set ray euler to the current ray euler (needed for sky shader)

            //run sky shader (wrapper)
            ray.color = sky_shader(ray); //sky color
            ray.reflection_color.empty=true; //set empty flag to true
            return ray;
        }

        // Select the closest intersect
        int closestIndex = -1;  // Initialize with an invalid index
        float closestDepth = std::numeric_limits<float>::max();  // Use the maximum possible depth initially

        for (int i = 0; i < intersects.size(); i++) {
            if (intersects[i].depth < closestDepth) {  // Check for the closest (smallest depth)
                closestDepth = intersects[i].depth;
                closestIndex = i;
            }
        }

        //assign the reference to the closest intersect to a variable
        Intersect& closestIntersect = intersects[closestIndex];


        //calculate frensel
        double frensel_deg = 180.0 - get_vectors_angle(closestIntersect.surface_normal, ray_euler);

        //TODO compute illumination...
        Lux lux = calculate_luminance(closestIntersect.intersect_point);

        //calculate reflection
        Vec3D reflection_vector = reflect_vector_normal(ray_euler, closestIntersect.surface_normal);

        //cast next bounce
        double epsilon = 1e-9; //advance ray forward by some epsilon
        Vec3D advanced_ray = closestIntersect.intersect_point + reflection_vector.set_magnitude(epsilon);
        
        //cast reflection ray from intersect
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
        

        //based on hit object, get shader inputs (inherit from object)
        ShaderInputs& shaderInputs = closestIntersect.target.shaderInputs;

        //run material shader
        std::vector<int> pixel = shader_wrapper(ray, lux, state.world, shaderInputs);

        //assign color to ray struct
        ray.color = pixel;

        //return
        return ray;
    };


};



struct Pin {
    Vec3D& target;
    Vec3D value1;
    Vec3D value2;
    int f1;
    int f2;

    int interpolatorFlag;

    // Constructor
    Pin(Vec3D& target, 
        const Vec3D& value1, 
        const Vec3D& value2, 
        int f1, int f2, 
        int interpolatorFlag)
        : target(target), value1(value1), value2(value2), f1(f1), f2(f2), interpolatorFlag(interpolatorFlag) {}
};


class Animator {
    private:
    Render& render;  // Store reference to Render

    int step = 0; //private step counter

    void updateAnimations() {
        //update render state based on Pins
        for (int i = 0; i < pins.size(); i++) {
            //update animation and assign to reference
            double t = (
                static_cast<double>(step) / 
                abs(static_cast<double>(pins[i].f2) - static_cast<double>(pins[i].f1)
            ) - 
            (
                static_cast<double>(pins[i].f1) / 
                abs(static_cast<double>(pins[i].f2) - static_cast<double>(pins[i].f1))
            ));
            
            //static_cast<double>(step) / static_cast<double>(pins[i].f2)
            Vec3D cache = vector_interpolate(
                pins[i].value1,
                pins[i].value2,
                t, //TODO correctly calculate
                pins[i].interpolatorFlag
            );

            //if pin is active apply animation
            if (pin_interpolator_range_checker(t, pins[i].interpolatorFlag)) {
                pins[i].target.x = cache.x;
                pins[i].target.y = cache.y;
                pins[i].target.z = cache.z;
            };

        };

    };

    public:
    // Constructor that takes a reference to Render
    Animator(Render& render) : render(render) {}

    //public variables
    int steps = 240;

    std::vector<Pin> pins={};

    //method declarations
    void addPin(Pin& pin) {
        pins.push_back(pin);
    };
    
    void RenderFrames(){
        for (int i = 0; i < steps; i++) {
            //apply animation Pins
            updateAnimations();

            //render frame
            render.render();

            //output buffers to file
            render.buffer.output_buffers_to_file("./tmp/f"+std::to_string(i)+".ppm", false);
            
            //complete frame
            step++;
            std::cout << "Rendered Frame: " << i << endl;
        }
    };
};

#endif