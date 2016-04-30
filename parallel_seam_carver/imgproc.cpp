/******************************************
** File: imgproc.cpp
** Authors: Augusto Blomer, Victor Konen
**
** Converts images to text lab for
** seam_carver and back
** See readme.txt
*******************************************/
#include <stdio.h>
#include <math.h>
#include <cfloat>
#include <string>

#include "cimg.h"
#include "sl_vector.h"

using namespace cimg_library;

void rgbToLab(const char * img_fn, const char * txt_fn) {
    FILE *f;
    int width;
    int height;

    //get the input image
    CImg<double> input(img_fn);
    width = input.width();
    height = input.height();

    //open the output file
    f = fopen(txt_fn, "w");
    if (f == NULL) {
        printf("Invalid output file.\n");
        return;
    }

    //write the dimensions
    fprintf(f, "%d\n%d\n", width, height);

    //write the depth and spectrum
    fprintf(f, "%d\n%d\n", input.depth(), input.spectrum());

    //convert input image to lab
    CImg<double> lab = input.RGBtoLab();

    //write pixel data to file
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            fprintf(f, "%.10f %.10f %.10f\n", lab(i, j, 0), lab(i, j, 1), lab(i, j, 2));
        }
    }

    //close the output file
    fclose(f);
}

void labToRgb(const char * img_fn, const char * txt_fn) {
    FILE *f;
    int width;
    int height;
    int depth;
    int spectrum;

    double templ, tempa, tempb;

    //open the input file
    f = fopen(txt_fn, "r");
    if (f == NULL) {
        printf("Invalid input file.\n");
        return;
    }

    //read in the dimensions
    fscanf(f, "%d\n%d\n", &width, &height);

    //read in the depth and spectrum
    fscanf(f, "%d\n%d\n", &depth, &spectrum);

    //read in the image info
    SlVector3 *image = new SlVector3[width * height];
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            fscanf(f, "%lf %lf %lf\n", &templ, &tempa, &tempb);
            image[i * height + j][0] = templ;
            image[i * height + j][1] = tempa;
            image[i * height + j][2] = tempb;
        }
    }

    //close the input file
    fclose(f);

    //store the output image
    CImg<double> output(width, height, depth, spectrum, 0);
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            output(i, j, 0) = image[i * height + j][0];
            output(i, j, 1) = image[i * height + j][1];
            output(i, j, 2) = image[i * height + j][2];
        }
    }

    //save the output image
    CImg<double> rgb = output.LabtoRGB();
    if (strstr(img_fn, "png")) {
        rgb.save_png(img_fn);        
    } else if (strstr(img_fn, "jpg")) {
        rgb.save_jpeg(img_fn);       
    }
}

int main(int argc, char *argv[]) {
    std::string process_type(argv[1]);
    bool type_switch;

    if (argc < 4) {
        printf("Usage: ./imgproc <option> <image_file> <text_file>\n");
        return 0;
    }

    if (process_type == "-lab") {
        type_switch = true;
    } else if (process_type == "-rgb") {
        type_switch = false;
    } else {
        printf("Error: Enter one option: -lab || -rgb\n");
        return 0;
    }

    if (type_switch) {
        rgbToLab(argv[2], argv[3]);
    } else {
        labToRgb(argv[2], argv[3]);
    }

    return 0;
}