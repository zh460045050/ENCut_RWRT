function [ucm2] = im2ucm(imgFile, outFile, times)

gPb_orient = globalPb(imgFile, outFile, times);
ucm2 = contours2ucm(gPb_orient, 'doubleSize');
save(outFile,'ucm2');
