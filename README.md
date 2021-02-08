# Online Rain/Snow Removal from Surveillance Videos (TIP2020)

This paper proposes a new online rain/snow removal method from surveillance videos by fully encoding the dynamic statistics of both rain/snow and background scenes in a video along time into the model, and realizing it with an online mode to make it potentially available to handle constantly coming streaming video sequence. 

[Minghan Li](https://scholar.google.com/citations?user=LhdBgMAAAAAJ&hl=en&oi=ao),
[Xiangyong Cao](https://scholar.google.com/citations?user=IePM9RsAAAAJ&hl=en),
[Qian Zhao](https://scholar.google.com/citations?user=vM6yGTEAAAAJ&hl=en),
[Lei Zhang](https://scholar.google.com/citations?user=tAK5l1IAAAAJ&hl=en&oi=ao),
[Deyu Meng](https://scholar.google.com/citations?user=an6w-64AAAAJ&hl=en&oi=ao)

Please go to the [Homepage](https://sites.google.com/view/onlinetmscsc/) to obtain more information about our work.

## Overview
The diagram of the proposed OTMS-CSC model implemented on a video with dynamic background. 

![fig](https://github.com/MinghanLi/OTMSCSC_matlab_2020/tree/master/figures/1_dynamic.png)

## Installation

1. Clone this repo

2. Install Matlab 

3. Compile GCO-v3.0
- Download [gco-v3.0 library](https://vision.cs.uwaterloo.ca/code/), and unzip the file.

- Start Matlab, and make gco-v3.0\matlab your working directory or add it to your path.

- To test your installation of GCO_MATLAB, run the gco-v3.0\matlab\GCO_UnitTest command.

4. Download dataset (NTURain or your own videos) and put the file into the input folder

5. Run demo.m

## Citation
Please cite our paper if you find anything helpful,

```
@article{Li2021OnlineRR,
  title={Online Rain/Snow Removal From Surveillance Videos},
  author={Minghan Li and Xiangyong Cao and Q. Zhao and L. Zhang and Deyu Meng},
  journal={IEEE Transactions on Image Processing},
  year={2021},
  volume={30},
  pages={2029-2044}
}
```

## Contact
For more information please contact liminghan0330@gmail.com

