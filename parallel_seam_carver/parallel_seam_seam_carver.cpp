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
    double *path_cost = (double *) malloc(initial_height * sizeof(double)); // path costs that will be sent
    int rows_per_process = current_height / numprocs;
    int start = rows_per_process * rank;
    int end = rows_per_process * (rank + 1);
    if (rank == numprocs - 1)
        end += (current_height % numprocs);
    
    //find the lowest cost seam by computing the lowest cost paths to each pixel
    for (int y = start; y < end; y++) {
        
        if (rank > 0 && y == start) {
            MPI_Recv(path_cost, current_height, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < current_width; i++)
                path_costs[i * initial_height + (y - 1)] = path_cost[i];
        }
        
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
        
        if (rank < (numprocs - 1) && y == (end - 1)) {
            for (int i = 0; i < current_width; i++)
                path_cost[i] = path_costs[i * initial_height + y];
            MPI_Send(path_cost, current_height, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
        }
    }
    
    int *recvcounts = new int[numprocs];
    int *displs = new int[numprocs];
    
    for (int i = 0; i < numprocs; i++)
        recvcounts[i] = rows_per_process * initial_height;
    recvcounts[numprocs - 1] += (current_height % numprocs) * initial_height;
    
    displs[0] = 0;
    for (int i = 1; i < numprocs; i++)
        displs[i] = displs[i - 1] + recvcounts[i - 1];
    
    int *previous = (int *) malloc(initial_width * initial_height * sizeof(int));
    
    for (int i = 0; i < current_height; i++) {
        for (int j = 0; j < current_width; j++) {
            previous[i * initial_height + j] = previous_x[j * initial_height + i];
        }
    }
    
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                       &previous[0], recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
    
    for (int i = 0; i < current_height; i++) {
        for (int j = 0; j < current_width; j++) {
            previous_x[i * initial_height + j] = previous[j * initial_height + i];
        }
    }
    
    //find the xcoord the lowest cost seam starts at the bottom of the current image
    int x_coord = 0;
    if (rank == numprocs - 1) {
        for (int x = 0; x < current_width; x++) {
            if (path_costs[x * initial_height + current_height - 1] < path_costs[x_coord * initial_height + current_height - 1]) {
                x_coord = x;
            }
        }
    }
    
    MPI_Bcast(&x_coord, 1, MPI_INT, numprocs - 1, MPI_COMM_WORLD);

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
    
    // split up work between processes
    double *path_cost = (double *) malloc(initial_height * sizeof(double)); // path costs that will be sent
    int columns_per_process = current_width / numprocs;
    int start = columns_per_process * rank;
    int end = columns_per_process * (rank + 1);
    if (rank == numprocs - 1)
        end += (current_width % numprocs);
    
    //find the lowest cost seam by computing the lowest cost paths to each pixel
    for (int x = start; x < end; x++) {
        
        if (rank > 0 && x == start) {
            MPI_Recv(path_cost, current_height, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < current_height; i++)
                path_costs[(x - 1) * initial_height + i] = path_cost[i];
        }
        
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
        
        if (rank < (numprocs - 1) && x == (end - 1)) {
            path_cost = &path_costs[x * initial_height];
            MPI_Send(path_cost, current_height, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
        }
    }
    
    int *recvcounts = new int[numprocs];
    int *displs = new int[numprocs];
    
    for (int i = 0; i < numprocs; i++)
        recvcounts[i] = columns_per_process * initial_height;
    recvcounts[numprocs - 1] += (current_width % numprocs) * initial_height;
    
    displs[0] = 0;
    for (int i = 1; i < numprocs; i++)
        displs[i] = displs[i - 1] + recvcounts[i - 1];
    
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                       &previous_y[0], recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
    
    //find the ycoord the lowest cost seam starts at the right of the current image
    int y_coord = 0;
    if (rank == numprocs - 1) {
        for (int y = 0; y < current_height; y++) {
            if (path_costs[(current_width - 1) * initial_height + y] < path_costs[(current_width - 1) * initial_height + y_coord]) {
                y_coord = y;
            }
        }
    }
    
    MPI_Bcast(&y_coord, 1, MPI_INT, numprocs - 1, MPI_COMM_WORLD);

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
