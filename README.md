# Piranhas
A Toolkit for creating Peripheral Architectures

This toolkit was used to create the peripheral architectures of the paper "Can Peripheral Representations Improve Clutter Metrics on Complex Scenes?" by Arturo Deza & Miguel Eckstein, 2016.

The toolkit was originally written in MATLAB, but has been extended to python and Torch to increase cross-collaborations between fields of vision science and computer vision and deep learning. With recent trends in RNN's, and Attention-based Networks, we decided to make our toolkit public to stimulate 'hybrid' ideas in the general vision (human,computer,robot) community.

# What is a Peripheral Architecture?
A peripheral architecture is a collection of regions that simulate human-like pooling regions and foveal and peripheral like mechanisms. Peripheral Architectures can be simple such as the ones used in "Multiple Object Recognition with Visual Attention" [Ba,Mnih & Kavukucoglu, 2015], or can be more complex given biological constraints such as the one proposed in "Metamers of the Ventral Stream" [Freeman & Simoncelli, 2011]. 

# Creating Peripheral Architectures

1) Download our toolbox to create Peripheral Architectures at your convenience in MATLAB, python or Torch.
2) Define your viewing parameters.

## Define your Computer parameters:

In python/piranhas.py @param_init():

	```python
  pixel_res_width = 800 # in pixels
  pixel_res_height = 600 # in pixels
  mon_width = 37.5 # in cm
  mon_height = 30 # in cm
	```

## Define your Human (Viewing Distance) parameters:

In python/piranhas.py @param_init():

	```python
  view_dist = 64 # in cm
  cm_per_pixel = mon_width/pixel_res_width;
  deg_per_pixel = 2*math.degrees(math.atan(cm_per_pixel/2/view_dist))
	```

