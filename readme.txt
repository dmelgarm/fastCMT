WARNING: This code is undergoing major overhaul to accomodate multiple point sources and is
changing frequently, use at your own risk


*****************************
********   fastCMT   ********
*****************************

VERSION 1.0

A program to compute moment tensors and locate centroids using coseismic offsets.

Diego Melgar 08/2011

Institue of Geophysics and Planetary Physics
Scripps Institution of Oceanograpy
University of California San Diego

dmelgarm@ucsd.edu
http://igpppublic.ucsd.edu/~dmelgarm

*****************************




You can do what you will with this code please acknowledge

* Melgar, D., Bock, Y. & Crowell, B.W., (2011), Real-Time Centroid Moment Tensor Determination for Large Earthquakes from Local and Regional Displacement Records, Geophys. J. Int., in review.

In that paper you can find more information on the ideas behind the code.


PROGRAM FLOW

Download the zip file and unpack it in your Matlab path you should see Matlab files

batch_fcmt
batch_green
buildrotmat
cmtgridmin
fastCMT
getmu
interleave
makegreen
makegrid
make_psmeca
mtinv2mt
prepare_elmayor
rgreen

and a GMT script

fastCMTtoki.gmt

Then proceed with the following steps

0. Read this readme.

1. Prepare the data. fastCMT uses data saved in a matlab structure with the following fields:

       N: North/South time series for all stations
       E: East/West time series for all stations
       U: Up/Down time series for all stations
       T: Time vector for all stations
     lat: Latitude of all stations
     lon: Longitude of all stations
    stdn: Standard deviation of pre-event noise in the north direction
    stde: Standard deviation of pre-event noise in the east direction
    stdu: Standard deviation of pre-event noise in the up direction

I have included a script prepare_elmayor() that I used to pre-process the El Mayor-Cucapah data, you probably can't use the same script for your data but you can modify it and it should serve you as a guide to what you need to do.

2. Get the EDGRN/EDCMP codes from: Wang, R., Martin, F.L. & Roth, F., 2003, Computation of deformation induced by earthquakes in a multi layered elastic crust – FORTRAN programs EDGRN/EDCMP, Computers & Geosciences, 29, 195-207. It can be downlaoded from the Computers and Geosciences web page. So, pick a velocity model for your problem and run EDGRN to make the Green functions. EDCMP is not used by fastCMT, but it's still a handy tool to have.

You should save the velocity model as a .mat file called velmod.mat in the same directory where you save the data. For example a 4 layer velocity model over a half space I used looks like this:

0	3500	1800	2200
1000	3500	1800	2200
1000	5000	2890	2650
11000	5000	2890	2650
11000	6500	3740	2870
23000	6500	3740	2870
23000	8100	4680	3300
1000000	0	0	0

Where column 1 is layer thickness, column 2 is vp, column 3 is vs and column 4 is rho.

3. fastCMT locates the centroid by running multiple ivnersions at different centroid locations, so you need to make GFs for each inversion node. Program batch_green() does this for you. You have to have run the EDGRN program for your velocity model and have GFs saved somewhere beforehand. The input variable 'data' is a structure with the fields listed in step 1 and it's made by a script like prepare_elmayor() as in step 1. This will make a .mat file with the GFs to be used in the next step.

4. Run the inversion. This is achieved with the program batch_fcmt() which invokes fastCMT() for each inversion node. The routine is commented and self-explanatory. This will output the inversion node coordinates, the moment tensor at each node at all epochs and the misfit at each node and all epochs. If you want to perform L1 inversions you need to download the L1 MAGIC suite of Matlab codes, Google it and you'll find it, it's free.

5. To analyze the data I have created a routine cmtgridmin() which picks out the best fitting moment tensor at each epoch and writes text files for plotting in GMT.

6. I have included an example of a GMT routine I use to make movies of the inversion results (for Tokachi-Oki in this case), it's called fastCMTtoki.gmt you should modify it to suit your needs. It uses ffmpeg to encode the movie, you can get that off the web.

Happy trails!

