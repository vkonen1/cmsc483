/******************************************
** File: seam_carver.cpp
** Authors: Augusto Blomer, Victor Konen
**
** Driver source code for the program
** See readme.txt
*******************************************/
#include <stdio.h>
#include <math.h>
#include <cfloat>
#include <string>

#include "cimg.h"
#include "sl_vector.h"

using namespace cimg_library;

//keeps track of the pixels in the image
SlVector3 *image;
int initial_width = 0;
int initial_height = 0;

int current_width = 0;
int current_height = 0;

int target_width = 0;
int target_height = 0;

//tracks the energy for each pixel
double *image_energy;
double max_energy = 0;

double *path_costs;
int *previous_x;
int *previous_y;

//computes the x gradient component for the energy function
double gradx(int x, int y) {
	double result = 0;
	if (x == 0) {
		result = 2.0 * mag(image[initial_height + y] - image[y]);
	} else if (x == initial_width - 1) {
		result = 2.0 * mag(image[x * initial_height + y] - image[(x - 1) * initial_height + y]);
	} else {
		result = mag(image[(x + 1) * initial_height + y] - image[(x - 1) * initial_height + y]);
	}

	if (result < 0) {
		result *= -1.0;
	}

	return result;
}

//computes the y gradient component for the energy function
double grady(int x, int y) {
	double result = 0;
	if (y == 0) {
		result = 2.0 * mag(image[x * initial_height + 1] - image[x * initial_height]);
	} else if (y == initial_height - 1) {
		result = 2.0 * mag(image[x * initial_height + y] - image[x * initial_height + y - 1]);
	} else {
		result = mag(image[x * initial_height + y + 1] - image[x * initial_height + y - 1]);
	}

	if (result < 0) {
		result *= -1.0;
	}

	return result;
}

//computes the energy of the pixel at x,y and updates the max energy if needed
double energy(int x, int y) {
	double result = sqrt(pow(gradx(x, y), 2) + pow(grady(x, y), 2));
	if (result > max_energy) {
		max_energy = result;
	}
	return result;
}

//finds the energy of each pixel in the image and stores it
void computeImageEnergy() {
	//compute the energy of the pixels
	for (int i = 0; i < initial_width; i++) {
		for (int j = 0; j < initial_height; j++) {
			image_energy[i * initial_height + j] = energy(i, j);
		}
	}
}

//removes the lowest energy vertical seam from the image
void removeVerticalSeam() {
	double energies[3];
	double min_energy;
	int prev_x;
	int prev_y;
	//find the lowest cost seam by computing the lowest cost paths to each pixel
	for (int y = 0; y < current_height; y++) {
		for (int x = 0; x < current_width; x++) {
			if (y == 0) {
				path_costs[x * initial_height] = image_energy[x * initial_height];
				previous_x[x * initial_height] = -1;
				previous_y[x * initial_height] = -1;
			} else {
				//the pixel directly above
				energies[1] = path_costs[x * initial_height + y - 1];
				//pixel above to the left
				if (x != 0) {
					energies[0] = path_costs[(x - 1) * initial_height + y - 1];
				} else {
					energies[0] = DBL_MAX;
				}
				//pixel above to the right
				if (x != current_width - 1) {
					energies[2] = path_costs[(x + 1) * initial_height + y - 1];
				} else {
					energies[2] = DBL_MAX;
				}

				//find the one with the least path cost
				min_energy = energies[0];
				prev_x = x - 1;
				prev_y = y - 1;
				if (energies[1] < min_energy) {
					min_energy = energies[1];
					prev_x = x;
				}
				if (energies[2] < min_energy) {
					min_energy = energies[2];
					prev_x = x + 1;
				}

				//set the minimum path cost for this pixel
				path_costs[x * initial_height + y] = min_energy + image_energy[x * initial_height + y];
				//set the previous pixel on the minimum path's coordinates for this pixel
				previous_x[x * initial_height + y] = prev_x;
				previous_y[x * initial_height + y] = prev_y;
			}
		}
	}

	//find the xcoord the lowest cost seam starts at the bottom of the current image
	int x_coord = 0;
	for (int x = 0; x < current_width; x++) {
		if (path_costs[x * initial_height + current_height - 1] < path_costs[x_coord * initial_height + current_height - 1]) {
			x_coord = x;
		}
	}

	//delete the seam from the bottom up
	for (int y = current_height - 1; y >= 0; y--) {
		//delete this pixel by copying over it and all those following to the right
		for (int x = x_coord; x < current_width - 1; x++) {
			image[x * initial_height + y] = image[(x + 1) * initial_height + y];
		}
		//next pixel
		x_coord = previous_x[x_coord * initial_height + y];
	}

	//decrease the current width of the image
	current_width--;
}

//removes the lowest energy vertical seam from the image
void removeHorizontalSeam() {
	double energies[3];
	double min_energy;
	int prev_x;
	int prev_y;
	//find the lowest cost seam by computing the lowest cost paths to each pixel
	for (int x = 0; x < current_width; x++) {
		for (int y = 0; y < current_height; y++) {
			if (x == 0) {
				path_costs[x * initial_height + y] = image_energy[x * initial_height + y];
				previous_x[x * initial_height + y] = -1;
				previous_y[x * initial_height + y] = -1;
			} else {
				//the pixel directly left
				energies[1] = path_costs[(x - 1) * initial_height + y];
				//pixel left and above
				if (y != 0) {
					energies[0] = path_costs[(x - 1) * initial_height + y - 1];
				} else {
					energies[0] = DBL_MAX;
				}
				//pixel left and below
				if (y != current_height - 1) {
					energies[2] = path_costs[(x - 1) * initial_height + y + 1];
				} else {
					energies[2] = DBL_MAX;
				}

				//find the one with the least path cost
				min_energy = energies[0];
				prev_x = x - 1;
				prev_y = y - 1;
				if (energies[1] < min_energy) {
					min_energy = energies[1];
					prev_y = y;
				}
				if (energies[2] < min_energy) {
					min_energy = energies[2];
					prev_y = y + 1;
				}

				//set the minimum path cost for this pixel
				path_costs[x * initial_height + y] = min_energy + image_energy[x * initial_height + y];
				//set the previous pixel on the minimum path's coordinates for this pixel
				previous_x[x * initial_height + y] = prev_x;
				previous_y[x * initial_height + y] = prev_y;
			}
		}
	}

	//find the ycoord the lowest cost seam starts at the right of the current image
	int y_coord = 0;
	for (int y = 0; y < current_height; y++) {
		if (path_costs[(current_width - 1) * initial_height + y] < path_costs[(current_width - 1) * initial_height + y_coord]) {
			y_coord = y;
		}
	}

	//delete the seam from right to left
	for (int x = current_width - 1; x >= 0; x--) {
		//delete this pixel by copying over it and all those following to the bottom
		for (int y = y_coord; y < current_height - 1; y++) {
			image[x * initial_height + y] = image[x * initial_height + y + 1];
		}
		//next pixel
		y_coord = previous_y[x * initial_height + y_coord];
	}

	//decrease the current height of the image
	current_height--;
}

//outputs a greyscale version of the input image based on the energy function
//used for testing
void outputGreyscale(char *file_name, double depth, double spectrum) {
	CImg<double> output(initial_width, initial_height, depth, spectrum, 0);
	for (int i = 0; i < initial_width; i++) {
		for (int j = 0; j < initial_height; j++) {
			double value = pow(image_energy[i * initial_height + j] / max_energy, 1.0 / 3.0) * 255.0;
			output(i, j, 0) = value;
			output(i, j, 1) = value;
			output(i, j, 2) = value;
		}
	}

	//save the output image
	CImg<double> rgb = output;
	if (strstr(file_name, "png")) {
		rgb.save_png(file_name);		
	} else if (strstr(file_name, "jpg")) {
		rgb.save_jpeg(file_name);		
	}
}

//outputs the seam carved version of the image at the new dimensions
void outputCarved(char *file_name, double depth, double spectrum) {
	//store the output image
	CImg<double> output(target_width, target_height, depth, spectrum, 0);
	for (int i = 0; i < target_width; i++) {
		for (int j = 0; j < target_height; j++) {
			output(i, j, 0) = image[i * initial_height + j][0];
			output(i, j, 1) = image[i * initial_height + j][1];
			output(i, j, 2) = image[i * initial_height + j][2];
		}
	}

	//save the output image
	CImg<double> rgb = output.LabtoRGB();
	if (strstr(file_name, "png")) {
		rgb.save_png(file_name);		
	} else if (strstr(file_name, "jpg")) {
		rgb.save_jpeg(file_name);		
	}
}

//stores the input image, resizes with seam carving, and generates the output image
int main(int argc, char *argv[]) {
	std::string process_type(argv[1]);
	bool carve;

	if (process_type == "-carve") {
		carve = true;
	} else if (process_type == "-greyscale") {
		carve = false;
	} else {
		printf("Error: Enter one option: -carve | -greyscale\n");
		return 0;
	}

	//get the input image
	CImg<double> input(argv[2]);
	initial_width = input.width();
	initial_height = input.height();

	current_width = initial_width;
	current_height = initial_height;

	//get the output width and height
	target_width = atoi(argv[4]);
	target_height = atoi(argv[5]);

	//exit if the output dimensions are larger than the input image dimensions
	if (initial_width < target_width || initial_height < target_height) {
		printf("Error: Desired output image dimensions larger than input image dimensions.\n");
		return 0;
	}

	//store the input image
	CImg<double> lab = input.RGBtoLab();
	image = new SlVector3[initial_width * initial_height];
	for (int i = 0; i < initial_width; i++) {
		for (int j = 0; j < initial_height; j++) {
			image[i * initial_height + j][0] = lab(i, j, 0);
			image[i * initial_height + j][1] = lab(i, j, 1);
			image[i * initial_height + j][2] = lab(i, j, 2);
		}
	}

	image_energy = new double[initial_width * initial_height];
	path_costs = new double[initial_width * initial_height];
	previous_x = new int[initial_width * initial_height];
	previous_y = new int[initial_width * initial_height];

	if (carve) {
		//remove vertical seams until reaching the target width
		while (current_width > target_width) {
			//calculate the energy of the image
			computeImageEnergy();
			//remove the seam
			removeVerticalSeam();
		}

		//remove horizontal seams until reaching the target height
		while (current_height > target_height) {
			//calculate the energy of the image
			computeImageEnergy();
			//remove the seam
			removeHorizontalSeam();
		}

		outputCarved(argv[3], input.depth(), input.spectrum());		
	} else {
		computeImageEnergy();
		outputGreyscale(argv[3], input.depth(), input.spectrum());
	}


	//clean up
	delete [] image;
	delete [] image_energy;
	return 0;
}
