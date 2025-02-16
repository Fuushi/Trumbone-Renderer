#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <math.h>

using namespace std;

#include "lib.h"


//main function
int main() {

    //create mesh generator
    Procedural_meshes generator;
    
    //create state object
    State state;

    //configure camera
    //state.camera.vec3D_c_vec = {0,-1.0,0};
    //state.camera.vec3D_c_pos = {3, 10, 0};
    
    //configure render options
    state.camera.res={512,512};
    state.camera.max_bounce=0;

    //initialize buffer
    OutputBuffer buffer;
    buffer.initializeBuffers(state.camera.res);

    //create renderer
    Render camera(state, buffer);

    //create first element, (plane)
    Mesh myMesh;
    generator.plane(myMesh);
    Element element(std::move(myMesh));
    //element.shaderInputs.material_color={255,0,0};
    element.shaderInputs.objectID=1;
    state.world.addElement(element);

    //create element 2 (cube)
    Mesh cubeMesh;
    generator.cube(cubeMesh, 5.0);
    Element element2(std::move(cubeMesh));
    //element2.shaderInputs.material_color={0,0,255};
    element.shaderInputs.objectID=2;
    state.world.addElement(element2);





    //create debug elements
    Mesh cubeMesh2;
    generator.cube(cubeMesh2, 0.5);
    Element element_x(std::move(cubeMesh2));
    element_x.shaderInputs.objectID=3;
    element_x.shaderInputs.material_color={255,0,0};
    element_x.pos={4,0,0};
    state.world.addElement(element_x);

    Mesh cubeMesh3;
    generator.cube(cubeMesh3, 0.5);
    Element element_y(std::move(cubeMesh3));
    element_y.shaderInputs.objectID=4;
    element_y.shaderInputs.material_color={0,255,0};
    element_y.pos={0,4,0};
    state.world.addElement(element_y);

    Mesh cubeMesh4;
    generator.cube(cubeMesh4, 0.5);
    Element element_z(std::move(cubeMesh4));
    element_z.shaderInputs.objectID=5;
    element_z.shaderInputs.material_color={0,0,255};
    element_z.pos={0,0,4};
    state.world.addElement(element_z);





    //create lights
    //L1
    Light light; //(sun)
    light.brightness=1;
    light.vec={-1.0,-1.0,-1.0};
    light.color={255,200,200};
    state.world.addLight(light);

    //L2
    Light light2;               //(Madara: "Oh? You Calculated Luminance for my light source?")
    light2.brightness=1;        //(Madara: "What about the second one?")
    light2.vec={-1.0,-1.0,-2.0};
    light2.color={200,200,255};
    state.world.addLight(light2);

    Animator animator(camera);
    animator.steps=60;

    //create Pin for camera pos
    Pin c_pin(
        state.camera.vec3D_c_vec,
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        0,
        60,
        0
    );

    Pin p_pin(
        state.camera.vec3D_c_pos,
        {3.0, -10.0, 0.0},
        {3.0, 0, -10.0},
        0,
        60,
        0
    );

    animator.addPin(c_pin);
    animator.addPin(p_pin);

    animator.RenderFrames();


    //begin render process

    //generate UVs
    //camera.render();
    
    //dump to file
    //buffer.output_buffers_to_file("./tmp/");

    //animator
    /*
    for (double i = 0; i < 60*10; i=i+1.0) {
                //generate UVs
        //camera.generate_uvs();
        //std::vector<double> rotation_uv_lmfao = {0.5, 0.6};
        //state.camera.vec3D_c_vec = rotation_matrix_degrees(
        //    state.camera.vec3D_c_vec, 75.0, rotation_uv_lmfao
        //);


        //camera.assign_rays();
        state.camera.vec3D_c_pos[0] = state.camera.vec3D_c_pos[0]+0.1;
        camera.render();
        //dump to file
        buffer.output_buffers_to_file("./tmp/");
    }
    */
    cout << "exit 1" << flush;
    return 1;
    
}