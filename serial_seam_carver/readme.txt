Authors: Augusto Blomer, Victor Konen

Code Notes:
Library Files:
CImg.H (Source included in header comment of source code) http://www.cimg.eu/
	Used to process images
	Convert images to raw data and back
slVector.H; slVector.cpp; slIO.H (Source included in header comments of source code)
	Used to store pixel color data
	Could be replaced with something simpler

Program Files:
seam_carver.cpp

Compilation Dependencies:
xorg-dev
libx11-dev
imagemagick
ghostscript

Overview:
Program takes in an input image and outputs an image that has been shrinked to the given dimensions using seam carving

Information on seam carving can be found here:
https://en.wikipedia.org/wiki/Seam_carving

Notes:
Assumes input image and output image are valid files
Exits gracefully if output_width or output_height are greater than the input image's corresponding dimensions

Building and Running:
Use 'make' to compile
Use './seam_carver option input_image output_image output_width output_height' to run

option is one of:
-carve
-greyscale