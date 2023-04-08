# Raytracing Algorithm
Ray-tracing based renderer from scratch for Computer Graphics course at DHBW Mannheim

There are a lot of different files here, so for better understanding:

* Documentation is stored as .docx Word-File and .pdf-File

* algo.cpp is the main file which is compiled and contains the main()-function

* because of poor performance in visual studio, the .exe files are run with cmd and the output of the program is stored to a ppm file with the `>` operator

* taskTwo.h is a headerfile which contains all the classes and different structs and methods which algo.cpp uses. *It is referenced by the path of the file in algo.cpp, so that probably has to be adjusted if this code is to be run locally.*

* taskOne.cpp is the main file which was used for the first task. After the second task, I quit using different files for different tasks for continuity.

Every image file is stored in .png and .ppm, while ppm was the output of the algorithm. PNG was needed to review pictures and import them to documentation.

* TaskOne is result image of task one, while TaskTwo is result of task two.

* sphere_one, two and five are self explanatory. These pictures contain a certain number of spheres.

* sphere_normal_shading is spheres rendered with color calculation by using normal vectors

* sphere_AAx1 to 16 are five spheres with Anti-Aliasing from 1-16 rays per pixel

* sphere_Diffuse_1 to 16 are spheres with diffuse material with recursion limit 1-16

* sphere_Diffuse_64 is a test with recursion limit 64 and different albedo

* sphere_Specular_1 to 16 are 5 spheres with diffuse and 3 with specular material with recursion limit 1-16

* sphere_Extrene are render tests with extremely high values for Anti-Aliasing and recursion limit, while the HD version is 1280x720 resolution

