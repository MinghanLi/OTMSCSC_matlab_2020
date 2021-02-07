# Online Rain/Snow Removal from Surveillance Videos

This paper proposes a new online rain/snow removal method from surveillance videos by fully encoding the dynamic statistics of both rain/snow and background scenes in a video along time into the model, and realizing it with an online mode to make it potentially available to handle constantly coming streaming video sequence. Specifically, inspired by the multi-scale convolutional sparse coding (MS-CSC) model designed for video rain removal (still for static rain), which finely delivers the sparse scattering and multi-scale shapes of real rain, this work encodes the dynamic temporal changing tendency of rain/snow and background motions as a dynamic MS-CSC framework by timely parameter amelioration for the model in an online implementation manner. Besides, a transformation operator capable of being adaptively updated along time is imposed on the background scenes to finely fit the background transformations existed in a video sequence. All these knowledge are formulated into a concise maximum a posterior (MAP) framework, which can be easily solved by alternative optimization technique.

## Overview
The diagram of the proposed OTMS-CSC model implemented on a video with dynamic background. As shown in the left figure, the background alignment based on its adjacent frames produces an initial stationary background. The online MS-CSC model shown on the right decomposes (a) the input video frame into four parts: (b) stationary background, (c) rain layer, (d) moving objects and the background noise. The rain layer (b) can be further decomposed as four sub-layers with various filters, which encode the repetitive local patterns of both rain/snow and background motions, displayed in the top-left corner of the second row. For the video with dynamic background, the final dynamic background (e) is the combination of the stationary background and the sub-layers with background motions, and the rectified rain layer only combines those sub-layers with relatively vertical filters representing rains.

![](https://github.com/MinghanLi/OTMSCSC_matlab_2020/tree/master/figures/1_dynamic.png)

## Requirements 
1. Matlab
2. GCO-v3.0
 
