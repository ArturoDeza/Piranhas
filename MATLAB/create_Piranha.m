%Set Initial Parameters:

function foveal_param = create_Piranha

% We Create a Peripheral Architecture with pooling regions 
% that approximate V1 receptive field size.

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

cm_per_pixel = param.mon_width/param.pixel_res_width;
deg_per_pixel = 2*atand(cm_per_pixel/2/param.view_dist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gabor Filter Bank Parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

orien = 8;
wave_num = 3;
%freq_zero = 2.0;
%freq_zero = 0.5;
freq_zero = 1.0;
bandwidth = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field of View Model Paramaters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Values used in the Deza & Eckstein, 2016 paper 
% Warning takes a bunch of space of memory!
%visual_field_radius_in_deg = 24;
%fovea = 2.0;
%scale = 0.25; %

visual_field_radius_in_deg = 14;
fovea = 1.0;
scale = 0.25;

e0_in_deg = 0.25;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization Parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

visual = 1;
visual_mask = 1;

% Compute other parameters:
[N_e N_theta] = get_pooling_parameters(scale,e0_in_deg,visual_field_radius_in_deg,deg_per_pixel);

% You can also manually put these in:
% N_e = 5;
% N_theta = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Peripheral Filters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

peripheral_filters = generate_pooling_regions_vector_smooth(deg_per_pixel, N_e, N_theta, visual_field_radius_in_deg, fovea,e0_in_deg, visual);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Peripheral Mask stream %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Nt_check Ne_check peri_height peri_width] = size(peripheral_filters.regions);

N_theta = Nt_check;
N_e = Ne_check;

select_mask_stream = cell(N_theta*N_e,1);
select_mask_label = zeros(peri_height,peri_width);

for z = 1:N_theta*N_e
	i1 = rem(z,N_theta);
	if i1==0
		i1 = N_theta;
	end

	j1 = floor((z-0.1)/N_theta)+1;
	select_mask_stream{z} = squeeze(peripheral_filters.regions(i1,j1,:,:));

	%New implementation:
	select_mask_label(sub2ind(size(select_mask_stream{z}),find(select_mask_stream{z})))=z;
end

%%%%%%%%%%%%%%%%%%%
% Get Foveal Mask %
%%%%%%%%%%%%%%%%%%%

for i=1:N_theta
	for j=1:size(peripheral_filters.offsets{i,1})
		point_buff = peripheral_filters.offsets{i,1}(j,:);
		point_vector(j) = norm(point_buff);
	end
	dist_mask(i) = min(point_vector);
	clear point_vector;
end

% Get Foveal Radius
foveal_radius = ceil(mean(dist_mask));

foveal_radius_px = foveal_radius;
foveal_radius_deg = foveal_radius * deg_per_pixel;

foveal_mask = zeros(peri_height,peri_width);

mid_point_w = round(peri_width/2);
mid_point_h = round(peri_height/2);


for i=1:peri_width
	for j=1:peri_height
		if sqrt((i-mid_point_w)^2+(j-mid_point_h)^2)<=foveal_radius
			foveal_mask(j,i) = 1;
		end
	end
end 


%%%%%%%%%%%%%%%%%%%%%%
% Create Filter Bank %
%%%%%%%%%%%%%%%%%%%%%%

filter_bank = create_gabor_bank_mod(orien,freq_zero,wave_num,bandwidth,param,visual);

for mm=1:wave_num*orien
	i1 = rem(mm,wave_num);
	if i1==0
		i1 = wave_num;
	end

	j1 = floor((mm-0.1)/wave_num)+1;

	filter_bank_stream{mm} = filter_bank{i1}{j1}{1};
end

for mm=1:wave_num*orien
	i1 = rem(mm,wave_num);
	if i1==0
		i1 = wave_num;
	end

	j1 = floor((mm-0.1)/wave_num)+1;

	filter_bank_stream{mm+wave_num*orien} = filter_bank{i1}{j1}{2};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create list of Foveal Parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

foveal_param.N_theta = N_theta;
foveal_param.N_e = N_e;
foveal_param.wave_num = wave_num;
foveal_param.orien = orien;
foveal_param.peri_height = peri_height;
foveal_param.peri_width = peri_width;

foveal_param.select_mask_stream = select_mask_stream;
foveal_param.select_mask_label = select_mask_label;
foveal_param.foveal_mask = foveal_mask;
foveal_param.foveal_radius_px = foveal_radius_px;
foveal_param.foveal_radius_deg = foveal_radius_deg;
foveal_param.filter_bank = filter_bank;
foveal_param.filter_bank_stream = filter_bank_stream;
foveal_param.peripheral_filters = peripheral_filters;

foveal_param.visual_mask = visual_mask;

foveal_param.deg_per_pixel = deg_per_pixel;

foveal_param.monitor.pixel_res_width = param.pixel_res_width;
foveal_param.monitor.pixel_res_height = param.pixel_res_height;
foveal_param.monitor.mon_width = param.mon_width;
foveal_param.monitor.mon_height = param.mon_height;
foveal_param.monitor.view_dist = param.view_dist;
