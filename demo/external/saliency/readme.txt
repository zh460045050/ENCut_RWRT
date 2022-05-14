This code is a beta version of Matlab implementation of the paper 'Saliency Detection Based on Integration of Boundary and Soft-segmentation'.

Thanks to J. van de Weijer, Th. Gevers, J-M Geusebroek "Boosting Color Saliency in Image Feature Detection"
(PAMI06), for providing the salient points detector.


Matlab script file to start:
  runsaliency.m: run the saliency measure

!!Note that:
(1) You have download the SLIC superpixel software to get the 'dat' file of your own picture. 
(2) The PB code is necessary to obtain the boundary image.
(3) Make sure that your system is 32 bit that the KDE algorithm can be used only on 32 bit system. 