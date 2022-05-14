# ENCut_RWRT
## Overview

The official code of the paper: [Explored Normalized Cut with Random Walk Refining Term for Image Segmentation][arxiv] (IEEE TIP)

## Getting Started

Our demo code provides in "demo" (corresponding source code are given in "source"):

1. the bipartition-based segmentation
2. the cluster-based segmentation
3. the contour-based hierarchical segmentation

To test the performance of our ENCut, Please run the demo.m:
1. the parameter "k" is the number of partition for the first two method.
2. the parameter "thr" is the threshold for the UCM (the third method). 
3. you can use any image to replace the parameter "img".


We also pre-generate the segmentation result for this test image and save it in the document "segresult":
1. the file with suffix '_cluster' is result of the cluster-based segmentation for the ENCut
2. the file with suffix '_fcluster' is the result of the cluster-based segmentation for the f-ENCut (ENCut with fast exploring method discussesd in Sec.4 in our paper)
3. the file with suffix '_bipart' is the result of the bipartition-based segmentation for the ENCut.
4. the file with suffix '_ucm' is the result of the contour-based hierarchical segmentation for the ENCut (with the gPB-OWT-UCM framework). 


## Citation

If you find our work useful in your research, please cite:

@article{ENCut,
  title={Explored Normalized Cut With Random Walk Refining Term for Image Segmentation},
  author={Zhu, Lei and Kang, Xuejing and Ye, Lizhu and Ming, Anlong},
  journal={IEEE Transactions on Image Processing},
  volume={31},
  pages={2893--2906},
  year={2022},
  publisher={IEEE}
}

[arxiv]: https://ieeexplore.ieee.org/abstract/document/9745758/
