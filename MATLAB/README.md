#Piranhas [MATLAB R2012a]

1. Download the Piranhas toolbox for MATLAB. Read the [Tutorial](https://github.com/ArturoDeza/Piranhas/tree/master/Tutorial) to learn more about these parameters.
2. Define your Computer + Human perception parameters.
	```matlab
	% All of this code is implemented inside the create_Piranha function
	% showed in Step 3.
		
	%%%%%%%%%%%%%%%%%%%%%%
	% Monitor Parameters %
	%%%%%%%%%%%%%%%%%%%%%%

	param.pixel_res_width = 800; % in pixels
	param.pixel_res_height = 600; % in pixels

	%param.pixel_res_width = 1024;
	%param.pixel_res_height = 768;

	%param.pixel_res_width = 1280;
	%param.pixel_res_height = 1024;

	param.mon_width = 37.5; % in cm
	param.mon_height = 30; % in cm

	%param.view_dist = 76; % in cm
	param.view_dist = 64;
	%param.view_dist = 38;
	param.gamma_c = 1;
	param.psi_c = 0;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Field of View Model Paramaters %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Values used in the Deza & Eckstein, 2016 paper 
	% Warning takes a bunch of space of memory!
	visual_field_radius_in_deg = 24;
	fovea = 2.0;
	scale = 0.25; %


	e0_in_deg = 0.25;

	visual_field_radius_in_deg = 10;
	fovea = 1.0;
	scale = 0.25;
	```


3. Create a Peripheral Architecture.
	```matlab
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Create Peripheral Architecture %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 
	foveal_param = create_Piranha;

	foveal_param_cheap.peri_height = foveal_param.peri_height;
	foveal_param_cheap.peri_width = foveal_param.peri_width;
	foveal_param_cheap.foveal_radius_px = foveal_param.foveal_radius_px;
	foveal_param_cheap.N_e = foveal_param.N_e;
	foveal_param_cheap.N_theta = foveal_param.N_theta;
	foveal_param_cheap.select_mask_stream = foveal_param.select_mask_stream;
	foveal_param_cheap.foveal_mask = foveal_param.foveal_mask;
	foveal_param_cheap.deg_per_pixel = foveal_param.deg_per_pixel;
	foveal_param_cheap.monitor = foveal_param.monitor;
	foveal_param_cheap.peripheral_filters = foveal_param.peripheral_filters;

	clear foveal_param;

	% Get Peripheral and Foveal masks

	for i=1:length(foveal_param_cheap.select_mask_stream)
		peri_indx{i} = find(foveal_param_cheap.select_mask_stream{i}~=0);
	end

	foveal_indx = find(foveal_param_cheap.foveal_mask~=0);
	```

	If you run:
	```matlab
	>>> piranha_sample = create_Piranha_Bench;
	```
	
	you should see the following:

	![PirArchMATLAB](http://imgur.com/SkHypjR.png)

4. Pool your dense feature maps as in (Deza & Eckstein, 2016):

	```matlab

	% For this demo we are loading an image and then computing Feature Congestion.
	% You can directly load any type of dense map as well.

	img = imresize(imread('../images/1.jpg'),0.5);
	[fc_score dense_map] = getClutter_FC(img);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Run pool_dense_map function %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[size_map_h size_map_w] = size(dense_map);

	% Assume a center fixation for now, you can change this with any fixation point.

	y_coord = round(size_map_h/2);
	x_coord = round(size_map_w/2);

	visual_toggle = 0;

	% NaN padding?
	y_lower_in = [];
	y_upper_in = [];
	x_lower_in = [];
	x_upper_in = [];

	foveal_map = pool_dense_map(foveal_param_cheap,dense_map,peri_indx,foveal_indx,y_coord,x_coord,visual_toggle,y_lower_in,y_upper_in,x_lower_in,x_upper_in))
	```

	![DenseAndFoveatedFeatures](http://i.imgur.com/VL0g79x.png)

5. For a demo run from the MATLAB command line prompt:

	```matlab
	> pool_dense_map_demo1
	```

### Tips and tricks:
Save your peripheral architectures in a folder!
