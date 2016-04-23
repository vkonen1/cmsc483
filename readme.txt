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
img_proc.cpp

Compilation Dependencies:
xorg-dev
libx11-dev
imagemagick
ghostscript

Overview:
Programs take in an input image and outputs an image that has been shrinked to the given dimensions using seam carving

Information on seam carving can be found here:
https://en.wikipedia.org/wiki/Seam_carving

Notes:
Assumes input image and output image are valid files
Exits gracefully if output_width or output_height are greater than the input image's corresponding dimensions

Building and Running:
Use 'make' to compile programs

Use './imgproc option image_file.(png|jpg) text_file.txt' to convert images to lab or back
option is -rgb to convert lab text file to image
option is -lab to convert image to lab text file

Use './seam_carver input_image_lab.txt output_image_lab.txt output_width output_height' to run