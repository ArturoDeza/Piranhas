#Piranhas [Torch7]

1. Download the Piranhas toolbox for torch. Read the [Tutorial](https://github.com/ArturoDeza/Piranhas/tree/master/Tutorial) to learn about these parameters.

2. Define your Computer + Human perception parameters.
	```lua
	function param:param_init_all(o)
	o = o or {}
	setmetatable(o, self)
	self.__index = self
	-- Computer parameters:
	self.pixel_res_width = 800 -- in pixels
	self.pixel_res_height = 800 -- in pixels
	self.mon_width = 37.5 -- in cm
	self.mon_height = 30 -- in cm
	self.view_dist = 64 -- in cm
	-- Computed from above:
	self.cm_per_pixel = (self.mon_width)/(self.pixel_res_width)
	self.deg_per_pixel = 2*math.deg(math.atan(self.cm_per_pixel/2/self.view_dist))
	-- Gabor Filter parameters
	self.orien = 8
	self.wave_num = 3
	self.freq_zero = 1.0
	self.bandwidth = 1.0
	self.gamma_c = 1 -- gabor parameter
	self.psi_c = 0 -- gabor parameter
	-- Human parameters:
	-- Field of View Model Parameters:
	self.visual_field_radius_in_deg = 10.0
	self.fovea = 1.0
	self.scale = 0.5
	-- e_0 set to 0.25
	self.e0_in_deg = 0.25
	-- visualization toggles
	self.visual = 1
	self.visual_mask = 1
	return o
	end
	```

3. Create a Peripheral Architecture
	```lua
	-- Create_Pianhas.lua script:

	-- Load the piranhas module
	require 'piranhas'
	-- require 'torch' (used inside piranhas.lua functions)

	param = param:param_init_all(nil)

	scale = param.scale
	fovea = param.fovea
	e0_in_deg = param.e0_in_deg
	visual_field_radius_in_deg = param.visual_field_radius_in_deg
	deg_per_pixel = param.deg_per_pixel

	-- Get Peripheral Architecture Parameters:

	N_e, N_theta  = get_pooling_parameters(scale,e0_in_deg,visual_field_radius_in_deg,deg_per_pixel)

	e_max = visual_field_radius_in_deg
	visual_field_width = math.floor(0.5 + 2*(visual_field_radius_in_deg/deg_per_pixel))

	regions = create_regions_vector_smooth(e0_in_deg,e_max,visual_field_width,deg_per_pixel,N_theta,N_e)
	```

	You should see something like this:
	
	![PirArchTorch](http://imgur.com/gGmWcRJ.png)

	Foveal Trimming is currently being developed for the Torch version.

4. Pool your dense feature maps.

	Coming soon!

5. For a demo run the script in Step 3 or
	```lua
	th> create_Piranhas.lua


### Tips and tricks:
Save your peripheral architectures in a folder!
