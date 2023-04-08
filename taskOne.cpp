#include <iostream>


struct ppm { //Abstracts Pixels with three values into single "objects"
    int r;
    int g;
    int b;
};
 

const int image_width = 480;
const int image_height = 360;

ppm image_color[image_width][image_height]; //2D-Array of Pixels that is used for storing pixels

/**
 * Function that outputs the currently stored image
*/
void output_image() 
{
    std::cout << "P3 \n" << image_width << ' ' << image_height << "\n255\n";
    for (int j = image_height-1; j >= 0; --j) 
    {
        std::cerr << "\rPixel lines remaining: " << j << ' ' << std::flush; // progress indicator for user comfort
        for (int i = 0; i < image_width; ++i) 
        {
                std::cout << image_color[i][j].r << ' ' << image_color[i][j].g << ' ' << image_color[i][j].b << "\n";
        }
    }
}

/**
 * main function creates an image with constant color
 * and calls the output-function
*/
int main()
{
    for (int j = image_height-1; j >= 0; --j) 
    {
        for (int i = 0; i < image_width; ++i) 
        {
                image_color[i][j] = {45, 186, 181};
        }
    }
    output_image();
    return 0;
}

