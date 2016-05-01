Authors: Augusto Blomer, Victor Konen

Code Notes:
Library Files:
CImg.H (Source included in header comment of source code) http://www.cimg.eu/
    Used to process images
    Convert images to raw data and back
slVector.H; slVector.cpp; slIO.H (Source included in header comments of source
code)
    Used to store pixel color data
    Could be replaced with something simpler

Program Files:
parallel_energy_seam_carver.cpp
parallel_seam_seam_carver.cpp
seam_carver.cpp
img_proc.cpp

Compilation Dependencies for imgproc:
xorg-dev
libx11-dev
imagemagick
ghostscript

Overview:
Programs take in an input image and outputs an image that has been shrinked to
the given dimensions using seam carving

imgproc
Converts images to lab format output via text files and vice versa

*seam_carver
Uses lab text files from imgproc as image data and carves image to specified
dimensions

Information on seam carving can be found here:
https://en.wikipedia.org/wiki/Seam_carving

Notes:
Assumes input image and output image are valid files
Exits gracefully if output_width or output_height are greater than the input
image's corresponding dimensions

Building and Running:
Use 'make <program name>' to compile program components
of 'make' to compile the combined components

Use './imgproc option image_file.(png|jpg) text_file.txt' to convert images to
lab or back
option is -rgb to convert lab text file to image
option is -lab to convert image to lab text file

Use './seam_carver input_image_lab.txt output_image_lab.txt output_width 
output_height' to run