#Piranhas [Python 2.7]

1. Download the Piranhas toolbox for python.

2. Define your Computer + Human perception parameters.

	Computer perception:
	```python
	>>> import piranhas

	>>>pixel_res_width = 800 # in pixels
	>>>pixel_res_height = 600 # in pixels
	>>>mon_width = 37.5 # in cm
	>>>mon_height = 30 # in cm
	>>>view_dist = 64 # in cm

	computer_input = piranhas.monitor(pixel_res_width,pixel_res_height,mon_width,mon_height,view_dist)
	```

	Human perception:
	```python
	>>>visual_field_radius_in_deg = 10
	>>>fovea = 1.0
	>>>scale = 0.5

	>>>human_input = piranhas.param(monitor_input,visual_field_radius_in_deg,fovea,scale)
	```
3. Create a Peripheral Architecture.

	```python
	python -i create_Piranhas.py
	```
4. Pool your dense feature maps:

5. For a demo run from the MATLAB command line prompt:
