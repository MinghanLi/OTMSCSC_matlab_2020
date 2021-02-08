# OTMSCSC_matlab_2020
Online Rain/Snow Removal from Surveillance Videos

This paper proposes a new online rain/snow removal method from surveillance videos by fully encoding the dynamic statistics of both rain/snow and background scenes in a video along time into the model, and realizing it with an online mode to make it potentially available to handle constantly coming streaming video sequence. Specifically, inspired by the multi-scale convolutional sparse coding (MS-CSC) model designed for video rain removal (still for static rain), which finely delivers the sparse scattering and multi-scale shapes of real rain, this work encodes the dynamic temporal changing tendency of rain/snow and background motions as a dynamic MS-CSC framework by timely parameter amelioration for the model in an online implementation manner. Besides, a transformation operator capable of being adaptively updated along time is imposed on the background scenes to finely fit the background transformations existed in a video sequence. All these knowledge are formulated into a concise maximum a posterior (MAP) framework, which can be easily solved by alternative optimization technique.

Overview



1. GCO-v3.0

2. 
