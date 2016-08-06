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

These parameters are useful to put in for two main reasons: 
1. If you are performing EyeTracking experiments and/or psychophysics, then you can then
calculate equivalent image/target appearance in terms of d.v.a. (degrees of visual angle). 
2. If creating a Peripheral Architecture,
it is better to make sure that is is having some sort of biological plausibility. Example: If I want to define a human fovea to have a value ranging between 1-2 degrees (as constrainted biologically), then I would like that the fovea occupies a reasonable size in pixel space. This will make more sense when we cover the basics of degrees of visual angle.

## Define your Human (Viewing Distance) parameters:

In python/piranhas.py @param_init():

```python
view_dist = 64 # in cm
cm_per_pixel = mon_width/pixel_res_width;
deg_per_pixel = 2*math.degrees(math.atan(cm_per_pixel/2/view_dist))
```

In the @param_init() class, we define the `deg_per_pixel` which indicates a rate of how many degrees a human observer is seeing per pixel. Conversely, pixel_per_degree would be the opposite rate and seems a bit easier to intepret since we are using values above 0. Example:
In our study (Deza & Eckstein, 2016) the deg_per_pixel is 0.022 for our EyeTracking Experiments, but the deg_per_pixel conversion scale in 
our peripheral architecture is 0.042 (twice the size, and computed with the values scenes in the example above).

Increasing the `deg_per_pixel` rate in a simulation has its advantages, mainly the fact that we can now downsample or resize our image input image,
to preserve the total input space in the degrees space. It is easy to see that if an image is viewd at `X_deg` by `Y_deg` degrees of visual angle, then 

	X_deg = x_pixels * deg_per_pixel (on a monitor) and Y_deg = y_pixels * deg_per_pixel (on a monitor)

but equivalent ways to get the same `X_deg` by `Y_deg` dimensions are:

	X_deg = (x_pixels/2) * (deg_per_pixels*2) and Y_deg = (y_pixels/2) * (deg_per_pixels*2)

or

	X_deg = (x_pixels*2) * (deg_per_pixels/2) and Y_deg = (y_pixels*2) * (deg_per_pixels/2)

## Degrees of Visual Angle Illustrated:

![EyeTrackerSettings](/images/EyeTrackerSettings.jpg)

Drawing courtesy of [Katie Koelher](http://koehler.moonfruit.com/home/4580555573)

Here is an example of how degrees of visual angle are calculated, this is a standard practice in vision science when working with EyeTracking technology (See EyeLink 1000 SR Research, User Manual for more information).

## [FAQ]: Frequenty Asked Questions:

### Q: It seems like there are many possible settings to create a `deg_per_pixel rate`, which ones do I pick?

A: In general there will be some contraints that are inmutable such as the dimensions of the monitor you are using, such as the monitor width and height of the visible area. Another semi-flexible constraint can be viewing distance in the EyeTracker. It is a general norm to have the monitor between 50 to up to 80 cm of viewing distance. You can play with these parameters. Other free parameters include monitor resolution (however these are discrete values contingent on quality of monitors: 800x600, or 1280x1024 are classic settings). In radiology monitors of about
the same physical dimensions can go up to: 2096x2800. These are commonly called 5 Megapixel resolution monitors. As a good rule of thumb, creating
`deg_per_pixel` values that go from 0.022 up to 0.044 are pretty reasonble to simulate in a lab environment.

### Q: Should degrees of visual angle be a standard in Computer Vision?

A: No.

### Q: I do Deep Learning, what can I get out of Piranhas toolkit?

### Q: I do Computer Vision, what can I get out of Piranhas toolkit?

### Q: I do Vision Science, what can I get out of Piranhas toolkit?



