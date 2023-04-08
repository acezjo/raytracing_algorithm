#include <iostream>
#include <random>
#include "C:\Users\Johnny\Documents\Render Algo\taskTwo.h"




const int image_width = 1280;
const int image_height = 720;

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

float getRandomFloat() {
    // create static generator for random numbers
    static std::mt19937 generator(std::random_device{}());
    // distribute floating point numbers between 0 and 1
    std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
    return distribution(generator);
}

/**
 * main function renders the image with background color
 * and hit objects
*/
int main()
{
    specular* material = new specular(color(0.7,0.2,0.2));
    diffuse* mat_d = new diffuse(color(0.1, 0.2, 0.8));
    diffuse* mat_dd = new diffuse(color(0.2,0.2,0.2));
    emissive* mat_e = new emissive(color(0.1,0.1,0.5), 70.0);
    
    sphere sphere_list[12];
    sphere_list[0] = sphere(vector(0, 0, -4), 0.7, color(1, 0, 0), mat_d);
    sphere_list[1] = sphere(vector(0.5, -0.5, -1), 0.3, color(1, 0.7 , 0), mat_d);
    sphere_list[2] = sphere(vector(1.0, 0.5, -2), 0.5, color(1, 1, 0.5), mat_dd);
    sphere_list[3] = sphere(vector(-1.4, 0.2, -3), 0.4, color(0.5, 0.2, 0.8), mat_d);
    sphere_list[4] = sphere(vector(-0.3, -0.7, -3), 0.4, color(1, 0, 0.7), mat_d);
    sphere_list[5] = sphere(vector(-0.9, -0.7, -2), 0.5, color(0.1, 0.8, 1), mat_dd);
    sphere_list[6] = sphere(vector(0, 1, -3), 0.4, color(0.5, 0.5, 0.5), material);
    sphere_list[7] = sphere(vector(-1.2, 0.7, -2), 0.4, color(0.8, 0.3, 0.5), material);
    sphere_list[8] = sphere(vector(0, -0.5, -2), 0.35, color(0.2, 0.6, 0.4), material);
    sphere_list[9] = sphere(vector(0, -0.5, -5), 1, color(0.2, 0.6, 0.4), mat_e);
    sphere_list[10] = sphere(vector(-1, -0.5, -0.2), 0.7, color(0.2, 0.6, 0.4), mat_e);
    sphere_list[11] = sphere(vector(1, 1.2, -4), 0.7, color(0.2, 0.6, 0.4), mat_e);



    std::vector<sphere> spheres_vec(std::begin(sphere_list), std::end(sphere_list));

    camera Cam;
    const int AAPixel = 32;
    const float InvAAPixel = 1.0 / AAPixel;
    const int depth = 16000;
    

    for (int j = image_height-1; j >= 0; --j) 
    {
        for (int i = 0; i < image_width; ++i) 
        {
            for (int a = 0; a < AAPixel; ++a) {
                auto u = double(i) / (image_width - getRandomFloat());
                auto v = double(j) / (image_height - getRandomFloat()); 
                Ray ray = Cam.get_ray(u,v);
                image_color[i][j] += ray_color(ray, spheres_vec, depth);
            }
            image_color[i][j] *= InvAAPixel;
        }
    }
    output_image();
    return 0;
}

