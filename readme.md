# Trumbone Raytracing Renderer

A proof-of-concept C++ library for raytracing and rendering scenes with elements, lights, and a camera system. Designed as a portfolio project, this library demonstrates fundamental concepts of raytracing and scene management.

## Features

- **Scene Management**: Add and inspect elements (meshes) and lights in a 3D scene.
- **Camera System**: Generate UV coordinates, assign rays, and manage rendering.
- **Output Buffering**: Save rendered outputs to files for further processing or visualization.
- **Modular Design**: Components like `State`, `Render`, and `OutputBuffer` provide a flexible structure for extending functionality.

## Prerequisites

- **C++17 or later**: The library utilizes modern C++ features.
- **Compiler**: Ensure you have a C++ compiler installed (e.g., GCC, Clang, or MSVC).
- **Standard Libraries**: Includes STL headers such as `<iostream>`, `<vector>`, `<fstream>`, `<cmath>`, and `<iomanip>`.

## Getting Started

### Installation

1. Clone the repository:
   ```bash
   git clone <repository_url>
   cd <repository_folder>
   ```
2. Include the library header file (`lib.h`) in your project.

### Usage Example

Here’s a basic example to demonstrate how to set up a scene, add elements and lights, and render the output:

```cpp
#include <iostream>
#include "lib.h"

int main() {
    // Create mesh generator
    Procedural_meshes generator; // Placeholder for model loading

    // Initialize state object
    State state;

    // Initialize output buffer
    OutputBuffer buffer;
    buffer.initializeBuffers(state.camera.res);

    // Create renderer
    Render camera(state, buffer);

    // Add first element (plane)
    Mesh myMesh;
    generator.plane(myMesh); // Placeholder for actual model loading
    Element element(std::move(myMesh));
    state.world.addElement(element);

    // Add second element (cube)
    Mesh cubeMesh;
    generator.cube(cubeMesh, 5.0); // Placeholder for actual model loading
    Element element2(std::move(cubeMesh));
    state.world.addElement(element2);

    // Inspect the first element
    state.world.elements[0].inspect();

    // Add a light to the scene
    Light light;
    light.brightness=1;
    light.vec={-1.0,-1.0,-1.0};
    light.color={255,255,255};
    state.world.addLight(light);

    // Begin rendering
    camera.render();

    // Save the output to files
    buffer.output_buffers_to_file("./tmp/");

    std::cout << "Rendering complete." << std::endl;
    return 0;
}
```

### Output

Rendered files are saved in the specified output directory (e.g., `./tmp/`).

## Project Structure

- **`State`**: Manages the overall scene, including elements, lights, and camera.
- **`Render`**: Handles rendering operations such as UV generation and ray assignment.
- **`OutputBuffer`**: Manages output buffers and handles file saving.
- **`Procedural_meshes`**: A placeholder for generating procedural meshes (to be replaced by model loading in the future).

## Future Plans

- Replace procedural mesh generation with robust model loading.
- Add support for advanced lighting and shading models.
- Implement real-time rendering capabilities.

## Contributing

Contributions are welcome! If you’d like to contribute, please fork the repository, create a feature branch, and submit a pull request.

## License

[MIT License](LICENSE)

---

