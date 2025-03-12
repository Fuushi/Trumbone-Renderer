#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <math.h>

//functions.h
#include "functions.h"
#include "mathHelper.h"
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

