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
    
    //configure render options
    state.camera.max_bounce=1;

    //initialize buffer
    OutputBuffer buffer;
    buffer.initializeBuffers(state.camera.res);

    //create renderer
    Render camera(state, buffer);

    //create first element, (plane)
    Mesh myMesh;
    generator.plane(myMesh);
    Element element(std::move(myMesh));
    state.world.addElement(element);

    //create element 2 (cube)
    Mesh cubeMesh;
    generator.cube(cubeMesh, 5.0);
    Element element2(std::move(cubeMesh));
    state.world.addElement(element2);

    //create lights
    Light light; //(sun)
    light.brightness=128;
    light.vec={-1.0,-1.0,-1.0};
    //light.color={255,0,0};

    //add light to scene
    state.world.addLight(light);


    Light light2; //(Madara: "What about the second one?")
    light2.brightness=128;
    light2.vec={-1.0,-1.0,-2.0};
    //light2.color={0,0,255};

    //add light to scene
    state.world.addLight(light2);


    
    
    //begin render process

    //generate UVs
    camera.render();
    
    //dump to file
    buffer.output_buffers_to_file("./tmp/");

    return 1; //disables animation player


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