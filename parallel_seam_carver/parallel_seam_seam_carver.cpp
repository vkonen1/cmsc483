/******************************************
** File: parallel_seam_seam_carver.cpp
** Authors: Augusto Blomer, Victor Konen
**
** Parallelized lowest cost seam finder only
** See readme.txt
*******************************************/
#include "mpi.h"

#include <stdio.h>
#include <math.h>
#include <cfloat>
#include <string>
#include <sstream>

#include "sl_vector.h"

#define USE_MPI

//keeps track of the pixels in the image
SlVector3 *image;
int initial_width = 0;
int initial_height = 0;

int depth = 0;
int spectrum = 0;

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

//number of pixels in the original image
int n;

//mpi variables
int rank, numprocs;

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
    
    // split up work between processes
    double *my_path_costs;
    double *my_previous_x;
    double *my_previous_y;
    double *temp_path_costs;
    double *temp_previous_x;
    double *temp_previous_y;
    int my_cols = current_width / numprocs;
    int low_cols = my_cols;
    int extra_cols = current_width % numprocs;
    int start;
    int x_offset;
    int recv_cols;

    double left_end_cost, right_end_cost, temp_end_cost;

    if (rank < extra_cols) {
        my_cols++;
        start = rank * my_cols;
    } else {
        start = (extra_cols * (my_cols + 1)) + ((rank - extra_cols) * my_cols);
    }

    //printf("%d %d %d\n", rank, start, my_cols);
    
    my_path_costs = (double *) malloc(my_cols * current_height * sizeof(double));
    my_previous_x = (double *) malloc(my_cols * current_height * sizeof(double));
    my_previous_y = (double *) malloc(my_cols * current_height * sizeof(double));
    //find the lowest cost seam by computing the lowest cost paths to each pixel
    for (int y = 0; y < current_height; y++) {
        //compute the path costs for my columns     
        for (int x = start; x < start + my_cols; x++) {
            //printf("%d %d %d %d %d\n", rank, x, y, (x - start) * current_height + y, my_cols * current_height);
            if (y == 0) {
                path_costs[x * initial_height] = image_energy[x * initial_height];
                my_path_costs[(x - start) * current_height + y] = path_costs[x * initial_height];

                previous_x[x * initial_height] = -1;
                my_previous_x[(x - start) * current_height + y] = previous_x[x * initial_height];

                previous_y[x * initial_height] = -1;
                my_previous_y[(x - start) * current_height + y] = previous_y[x * initial_height];
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
                my_path_costs[(x - start) * current_height + y] = path_costs[x * initial_height + y];

                //set the previous pixel on the minimum path's coordinates for this pixel
                previous_x[x * initial_height + y] = prev_x;
                my_previous_x[(x - start) * current_height + y] = previous_x[x * initial_height + y];

                previous_y[x * initial_height + y] = prev_y;
                my_previous_y[(x - start) * current_height + y] = previous_y[x * initial_height + y];
            }
        }

        //send path cost needed to neighboring processes
        if (numprocs > 1) {
            if (rank != numprocs - 1) {
                //send rightmost cost to following process
                right_end_cost = path_costs[(start + my_cols - 1) * initial_height + y];        
                MPI_Send(&right_end_cost, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);

                //receive following process's leftmost cost
                MPI_Recv(&temp_end_cost, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                path_costs[(start + my_cols) * initial_height + y] = temp_end_cost;
            }
            if (rank != 0) {
                //send leftmost cost to preceding process
                left_end_cost = path_costs[start * initial_height + y];
                MPI_Send(&left_end_cost, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);

                //receive preceding process's rightmost cost
                MPI_Recv(&temp_end_cost, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                path_costs[(start - 1) * initial_height + y] = temp_end_cost;
            }            
        }
    }

    //update path costs and previous for all processes
    for (int i = 0; i < numprocs; i++) {
        if (rank == i) {
            continue;
        }

        if (i < extra_cols) {
            x_offset = i * (low_cols + 1);
            recv_cols = (low_cols + 1);
        } else {
            x_offset = (extra_cols * (low_cols + 1)) + ((i - extra_cols) * low_cols);
            recv_cols = low_cols;
        }

        //printf("%d %d\n", low_cols, extra_cols);
        //printf("%d %d %d\n", rank, x_offset, recv_cols);

        temp_path_costs = (double *) malloc(recv_cols * current_height * sizeof(double));
        temp_previous_x = (double *) malloc(recv_cols * current_height * sizeof(double));
        temp_previous_y = (double *) malloc(recv_cols * current_height * sizeof(double));
        MPI_Sendrecv(my_path_costs, my_cols * current_height, MPI_DOUBLE, i, 0, 
            temp_path_costs, recv_cols * current_height, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, 
            MPI_STATUS_IGNORE);
        MPI_Sendrecv(my_previous_x, my_cols * current_height, MPI_DOUBLE, i, 1, 
            temp_previous_x, recv_cols * current_height, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, 
            MPI_STATUS_IGNORE);
        MPI_Sendrecv(my_previous_y, my_cols * current_height, MPI_DOUBLE, i, 2, 
            temp_previous_y, recv_cols * current_height, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, 
            MPI_STATUS_IGNORE);
        /* problem is here */
        for (int j = 0; j < recv_cols * current_height; j++) {
            int x = x_offset + (j % recv_cols);
            int y = j / recv_cols;
            //printf("%d %d %d %d %d\n", rank, x, y, x * initial_height + y, recv_cols * current_height);
            //printf("%d\n", initial_height * initial_width);
            path_costs[x * initial_height + y] = temp_path_costs[(x - x_offset) * current_height + y];
            previous_x[x * initial_height + y] = temp_previous_x[(x - x_offset) * current_height + y];
            previous_y[x * initial_height + y] = temp_previous_y[(x - x_offset) * current_height + y];
        }
        free(temp_path_costs);
        free(temp_previous_x);
        free(temp_previous_y);
    }
    free(my_path_costs);
    free(my_previous_x);
    free(my_previous_y);

    //printf("here\n");
    
    //find the xcoord the lowest cost seam starts at the bottom of the current image
    int x_coord = 0;
    for (int x = 0; x < current_width; x++) {
        if (path_costs[x * initial_height + current_height - 1] < path_costs[x_coord * initial_height + current_height - 1]) {
            x_coord = x;
        }
    }

    //printf("here\n");


    //delete the seam from the bottom up
    for (int y = current_height - 1; y >= 0; y--) {
        //delete this pixel by copying over it and all those following to the right
        for (int x = x_coord; x < current_width - 1; x++) {
            image[x * initial_height + y] = image[(x + 1) * initial_height + y];
        }
        //next pixel
        //printf("%d\n", x_coord * initial_height + y);
        x_coord = previous_x[x_coord * initial_height + y];
        //printf("%d %d\n", rank, x_coord);
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
    
    // split up work between processes
    double *my_path_costs;
    double *my_previous_x;
    double *my_previous_y;
    double *temp_path_costs;
    double *temp_previous_x;
    double *temp_previous_y;
    int my_rows = current_height / numprocs;
    int low_rows = my_rows;
    int extra_rows = current_height % numprocs;
    int start;
    int y_offset;
    int recv_rows;

    double top_end_cost, bottom_end_cost, temp_end_cost;

    if (rank < extra_rows) {
        my_rows++;
        start = rank * my_rows;
    } else {
        start = (extra_rows * (my_rows + 1)) + ((rank - extra_rows) * my_rows);
    }
    
    my_path_costs = (double *) malloc(my_rows * current_width * sizeof(double));
    my_previous_x = (double *) malloc(my_rows * current_width * sizeof(double));
    my_previous_y = (double *) malloc(my_rows * current_width * sizeof(double));
    //find the lowest cost seam by computing the lowest cost paths to each pixel
    for (int x = 0; x < current_width; x++) {
        //compute the path costs for my rows
        for (int y = start; y < start + my_rows; y++) {

            if (x == 0) {
                path_costs[x * initial_height + y] = image_energy[x * initial_height + y];
                my_path_costs[(y - start) * current_width + x] = path_costs[x * initial_height + y];

                previous_x[x * initial_height + y] = -1;
                my_previous_x[(y - start) * current_width + x] = previous_x[x * initial_height + y];

                previous_y[x * initial_height + y] = -1;
                my_previous_y[(y - start) * current_width + x] = previous_y[x * initial_height + y];
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
                my_path_costs[(y - start) * current_width + x] = path_costs[x * initial_height + y];

                //set the previous pixel on the minimum path's coordinates for this pixel
                previous_x[x * initial_height + y] = prev_x;
                my_previous_x[(y - start) * current_width + x] = previous_x[x * initial_height + y];

                previous_y[x * initial_height + y] = prev_y;
                my_previous_y[(y - start) * current_width + x] = previous_y[x * initial_height + y];
            }
        }

        //send path cost needed to neighboring processes
        if (numprocs > 1) {
            if (rank != numprocs - 1) {
                //send bottom most cost to following process
                bottom_end_cost = path_costs[x * initial_height + (start + my_rows - 1)];        
                MPI_Send(&bottom_end_cost, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);

                //receive following process's top most cost
                MPI_Recv(&temp_end_cost, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                path_costs[x * initial_height + (start + my_rows)] = temp_end_cost;
            }
            if (rank != 0) {
                //send top most cost to preceding process
                top_end_cost = path_costs[x * initial_height + start];
                MPI_Send(&top_end_cost, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);

                //receive preceding process's bottom most cost
                MPI_Recv(&temp_end_cost, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                path_costs[x * initial_height + (start - 1)] = temp_end_cost;
            }            
        }
    }

    //update paths costs for all processes
    for (int i = 0; i < numprocs; i++) {
        if (rank == i) {
            continue;
        }

        if (i < extra_rows) {
            y_offset = i * (low_rows + 1);
            recv_rows = low_rows + 1;
        } else {
            y_offset = (extra_rows * (low_rows + 1)) + ((i - extra_rows) * low_rows);
            recv_rows = low_rows;
        }

        //printf("%d %d\n", low_rows, extra_rows);
        //printf("%d %d %d\n", rank, y_offset, recv_rows);

        temp_path_costs = (double *) malloc(recv_rows * current_width * sizeof(double));
        temp_previous_x = (double *) malloc(recv_rows * current_width * sizeof(double));
        temp_previous_y = (double *) malloc(recv_rows * current_width * sizeof(double));
        MPI_Sendrecv(my_path_costs, my_rows * current_width, MPI_DOUBLE, i, 0, 
            temp_path_costs, recv_rows * current_width, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, 
            MPI_STATUS_IGNORE);
        MPI_Sendrecv(my_previous_x, my_rows * current_width, MPI_DOUBLE, i, 1, 
            temp_previous_x, recv_rows * current_width, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, 
            MPI_STATUS_IGNORE);
        MPI_Sendrecv(my_previous_y, my_rows * current_width, MPI_DOUBLE, i, 2, 
            temp_previous_y, recv_rows * current_width, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, 
            MPI_STATUS_IGNORE);
    
        for (int j = 0; j < recv_rows; j++) {
            int x = j / recv_rows;
            int y = y_offset + (j % recv_rows);
            //printf("%d %d %d %d %d\n", rank, x, y, x * initial_height + y, recv_rows * current_width);
            //printf("%d\n", initial_height * initial_width);
            path_costs[x * initial_height + y] = temp_path_costs[(y - y_offset) * current_width + x];
            previous_x[x * initial_height + y] = temp_previous_x[(y - y_offset) * current_width + x];
            previous_y[x * initial_height + y] = temp_previous_y[(y - y_offset) * current_width + x];
        }
        free(temp_path_costs);
        free(temp_previous_x);
        free(temp_previous_y);
    }
    free(my_path_costs);
    free(my_previous_x);
    free(my_previous_y);

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

bool readInput(const char * file_name) {
    FILE *f;

    double templ, tempa, tempb;

    //open the input file
    f = fopen(file_name, "r");
    if (f == NULL) {
        printf("Invalid input file.\n");
        return false;
    }

    //read in the dimensions
    fscanf(f, "%d\n%d\n", &initial_width, &initial_height);
    current_width = initial_width;
    current_height = initial_height;

    //read in the depth and spectrum
    fscanf(f, "%d\n%d\n", &depth, &spectrum);

    //read in the image info
    image = new SlVector3[initial_width * initial_height];
    for (int i = 0; i < initial_width; i++) {
        for (int j = 0; j < initial_height; j++) {
            fscanf(f, "%lf %lf %lf\n", &templ, &tempa, &tempb);
            image[i * initial_height + j][0] = templ;
            image[i * initial_height + j][1] = tempa;
            image[i * initial_height + j][2] = tempb;
        }
    }

    //close the input file
    fclose(f);

    return true;
}

//outputs the seam carved version of the image at the new dimensions
void outputCarved(const char *file_name) {
    FILE *f;

    //open the output file
    f = fopen(file_name, "w");
    if (f == NULL) {
        printf("Invalid output file.\n");
        return;
    }

    //write the dimensions
    fprintf(f, "%d\n%d\n", target_width, target_height);

    //write the depth and spectrum
    fprintf(f, "%d\n%d\n", depth, spectrum);

    //write pixel data to file
    for (int i = 0; i < target_width; i++) {
        for (int j = 0; j < target_height; j++) {
            fprintf(f, "%.10f %.10f %.10f\n", image[i * initial_height + j][0],
                image[i * initial_height + j][1], image[i * initial_height + j][2]);
        }
    }

    //close the output file
    fclose(f);
}

//stores the input image, resizes with seam carving, and generates the output image
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!readInput(argv[1])) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //get the output width and height
    target_width = atoi(argv[3]);
    target_height = atoi(argv[4]);

    //exit if the output dimensions are larger than the input image dimensions
    if (initial_width < target_width || initial_height < target_height) {
        printf("Error: Desired output image dimensions larger than input image dimensions.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    image_energy = (double *) malloc(initial_width * initial_height * sizeof(double));
    path_costs = (double *) malloc(initial_width * initial_height * sizeof(double));
    previous_x = (int *) malloc(initial_width * initial_height * sizeof(int));
    previous_y = (int *) malloc(initial_width * initial_height * sizeof(int));
    n = initial_width * initial_height;

    if (image_energy == NULL || path_costs == NULL || previous_x == NULL || previous_y == NULL) {
        printf("problem");
    }

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

    if (rank == 0) {
        outputCarved(argv[2]);
    }

    free(image_energy);
    free(path_costs);
    free(previous_x);
    free(previous_y);

    MPI_Finalize();

    return 0;
}
