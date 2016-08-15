## Define your Computer parameters:

In python/piranhas.py @param_init():

```python
pixel_res_width = 800 # in pixels
pixel_res_height = 600 # in pixels
mon_width = 37.5 # in cm
mon_height = 30 # in cm
```

These parameters are useful for two main reasons: 

* If you are performing EyeTracking experiments and/or psychophysics, you can then
calculate equivalent image/target appearance in terms of d.v.a. (degrees of visual angle). 
* If creating a Peripheral Architecture, it is better to be sure that it has having some sort of biological plausibility. Example: If you want to define a human fovea to have a value ranging between 1-2 degrees (as constrained biologically), then you would like that the fovea occupies a reasonable size in pixel space. This will make more sense when we cover the basics of degrees of visual angle.

## Define your Human parameters:

In python/piranhas.py @param_init():

```python
view_dist = 64 # in cm
cm_per_pixel = mon_width/pixel_res_width;
deg_per_pixel = 2*math.degrees(math.atan(cm_per_pixel/2/view_dist))


  ##################################
  # Field of View Model Parameters #
  ##################################
  
  visual_field_radius_in_deg = 10 # This is the total visual field_radius in deg
  fovea = 1.0 # This is the approximate foveal radius.
  #scale = 0.25 # Use this for V1 receptive field size pooling regions.
  scale = 0.5 # Use this for V2 receptive field size pooling regions.
  # See Simoncelli & Freeman, 2011 (Supplementary Material) for more information on the scale parameter.

```

In the` @param_init()` class, we define the `deg_per_pixel` which indicates a rate of how many degrees of visual angle a human observer is seeing per pixel. Conversely, `pixel_per_degree` would be the opposite rate and seems a bit easier to intepret since we are using values above 0. Example:
In our study (Deza & Eckstein, 2016) the deg_per_pixel is 0.022 for our EyeTracking Experiments, but the `deg_per_pixel` conversion scale in 
our peripheral architecture is 0.042 (twice the size, and computed with the values scenes in the example above).

Increasing the `deg_per_pixel` rate in a simulation has its advantages, mainly the fact that we can now downsample or resize our image input image,
to preserve the input information in the degrees space. It is easy to see that if an image is viewed at `X_deg` by `Y_deg` degrees of visual angle, then 

	X_deg = x_pixels * deg_per_pixel (on a monitor) and Y_deg = y_pixels * deg_per_pixel (on a monitor)

but equivalent ways to get the same `X_deg` by `Y_deg` dimensions are:

	X_deg = (x_pixels/2) * (deg_per_pixels*2) and Y_deg = (y_pixels/2) * (deg_per_pixels*2)

or

	X_deg = (x_pixels*2) * (deg_per_pixels/2) and Y_deg = (y_pixels*2) * (deg_per_pixels/2)

### Degrees of Visual Angle Illustrated:

![EyeTrackerSettings](/images/EyeTrackerSettings.jpg)

Drawing courtesy of [Katie Koelher](http://koehler.moonfruit.com/home/4580555573)

Here is an example of how degrees of visual angle are calculated, this is a standard practice in vision science when working with EyeTracking technology (See [EyeLink 1000 SR Research, User Manual](http://sr-research.jp/support/EyeLink%201000%20User%20Manual%201.5.0.pdf) for more information).

Here is a diagram of the human retina and an illustration of how sensor density decreases as a function of eccentricity:

![Diagram](http://www.webexhibits.org/causesofcolor/images/content/26z.jpg)
Courtesy of  www.webexhibits.org

## Create a Peripheral Architecture

[Matlab](https://github.com/ArturoDeza/Piranhas/tree/master/MATLAB):
```matlab
>>> pirAr = create_Piranha
```

[Python](https://github.com/ArturoDeza/Piranhas/tree/master/python):
```
$ python create_Piranha.py
```

[Torch](https://github.com/ArturoDeza/Piranhas/tree/master/torch):
```
$ th create_Piranha.lua
-- Or:
> require 'piranhas' 
```

## Download Peripheral Architectures from Piranhas School:

Feature coming soon!

### Q: It seems like there are many possible settings to create a `deg_per_pixel` rate, which ones do I pick?

A: In general there will be some contraints that are inmutable such as the dimensions of the monitor you are using, such as the monitor width and height of the visible area (screen). Another semi-flexible constraint can be viewing distance in the EyeTracker. It is a general norm to have the monitor between 50 to up to 80 cm of viewing distance. You can play with these parameters. Other free parameters include monitor resolution (however these are discrete values contingent on quality of monitors: 800x600, or 1280x1024 are classic settings). In radiology monitors of about
the same physical dimensions can go up to: 2096x2800. These are commonly called 5 Megapixel resolution monitors, and are quite expensive. As a good rule of thumb, creating `deg_per_pixel` values that go from 0.022 up to 0.044 are pretty reasonable to simulate in a lab environment, given common
viewing distances and hardware monitor settings.

