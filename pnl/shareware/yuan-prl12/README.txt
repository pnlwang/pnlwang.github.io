This folder contains the Matlab code for the regression (factorization) based image segmentation algorithm described in the following paper:

"Image segmentation using local spectral histograms and linear regression," by J. Yuan, D. L. Wang, and R. Li, in Pattern Recognition Letters, vol. 33, pp. 615-622, 2012.

Here are the instructions on how to run it:


1) After unzip the file, change the Current Directory in Matlab to the unzipped folder "FactorizationSeg".

2) In the Matlab command prompt, type demoFactorizationSeg to see a demo.

3) You can try any of the functions

Main functions:

sugImg.m -- Given an image, produce an N1*N2*bn array containing 'bn' filter responses    

factorizationSeg.m -- Produce segmentation results in form of a label map

NNMFSeg.m -- Given the SVD based solution, produce segmentation results with nonnegativity constraint.

Other comments:

- This program segments gray level images.  For color or multispectral images, run function 'subImg' on each band and combine the arrays using Ig=cat(3,Ig1,Ig2,...). For better results, convert RGB color values to the L*a*b*color space.

- In the demo, 11 filters are used (intensity + 2 LoG + 8 Gabor). For images without complicated textures, only a few filters (e.g., intensity + 1 LoG) can work sufficiently well, which can further reduce computation time.  
