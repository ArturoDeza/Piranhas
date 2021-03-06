clc;close all;clear all;

addpath('./denseClutterFunctions/matlabPyrTools/');
addpath('./denseClutterFunctions/Feature_Congestion/');

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

for i=1:length(foveal_param_cheap.select_mask_stream)
	peri_indx{i} = find(foveal_param_cheap.select_mask_stream{i}~=0);
end

foveal_indx = find(foveal_param_cheap.foveal_mask~=0);

%%%%%%%%%%%%%%%%%%%
% Load Dense map: %
%%%%%%%%%%%%%%%%%%%

% Insert code/function to compute or load a dense map here.

%dense_map = load('../images/1.mat');
%dense_map = dense_map.ProtoTotal;

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

y_lower_in = [];
y_upper_in = [];
x_lower_in = [];
x_upper_in = [];

%time operations?
time_operation = 1;

if time_operation
	tic
	foveal_map = pool_dense_map(foveal_param_cheap,dense_map,peri_indx,foveal_indx,y_coord,x_coord,visual_toggle,y_lower_in,y_upper_in,x_lower_in,x_upper_in);
	toc
else
	foveal_map = pool_dense_map(foveal_param_cheap,dense_map,peri_indx,foveal_indx,y_coord,x_coord,visual_toggle,y_lower_in,y_upper_in,x_lower_in,x_upper_in);
end





