# Python Script to initialize a Peripheral Architecture:
# Developed by Arturo Deza as part of the Piranhas Toolkit
# Questions & bugs: deza @ dyns.ucsb.edu

import numpy as np
from piranhas import *
import math

# Run Intialization Parameters:
param = param_init_all()

scale = param.scale
fovea = param.fovea
e0_in_deg = param.e0_in_deg
visual_field_radius_in_deg = param.visual_field_radius_in_deg
deg_per_pixel = param.deg_per_pixel

# Get Peripheral Architecture Parameters

# Given the scale, compute the N_e and N_theta parameters.
# You can also uncomment this line and input any values for N_e and N_theta.
# Please refer to Freeman & Simoncelli, 2011 to check what values of scale
# are biologically plausible for V1 and V2 receptive field sizes.

N_param = get_pooling_parameters(scale,e0_in_deg,visual_field_radius_in_deg,deg_per_pixel)

N_e = N_param[0]
N_theta = N_param[1]

#visual flag?
visual = 0

# The peripheral_regions data structure allows you to do anything you want with the pooling regions.
peripheral_filters = generate_pooling_regios_vector_smooth(deg_per_pixel,N_e,N_theta,visual_field_radius_in_deg,fovea,e0_in_deg,visual)

##################################
# Compute Peripheral Mask Stream #
##################################


# size might have changed after removing peripheral regions inside fovea
param_check = np.shape(peripheral_filters.regions)

N_theta = param_check[0]
N_e = param_check[1]
peri_height = param_check[2]
peri_width = param_check[3]

select_mask_stream = np.empty((N_theta*N_e,peri_height,peri_width),dtype=object)
select_mask_label = np.zeros((peri_height,peri_width))

for z in range(1,N_theta*N_e):
	i1 = np.mod(z,N_theta)
	if i1==0:
		i1 = N_theta
	
	j1 = np.round(np.floor((z-0.1)/N_theta)+1)
	select_mask_stream[z-1,:,:] = np.squeeze(peripheral_filters.regions[i1-1,j1-1,:,:])
	
	select_loc = np.where(select_mask_stream[z-1,:,:]>0)
	select_mask_label[select_loc[0],select_loc[1]] = z-1

###################
# Get Foveal Mask #
###################

dist_mask = np.empty(N_theta)
for i in range(0,N_theta):
	point_vector = np.empty(np.shape(peripheral_filters.offsets_r[i][0])[0])
	for j in range(0,np.shape(peripheral_filters.offsets_r[i][0])[0]):
		#Compute L2 norm (euclidean distance from center of fixation)
		point_buff = np.linalg.norm([peripheral_filters.offsets_r[i][0][j],peripheral_filters.offsets_c[i][0][j]])
		point_vector[j] = abs(point_buff)
	dist_mask[i] = min(point_vector)
	del point_vector

# Compute Foveal Radius
foveal_radius = np.ceil(np.mean(dist_mask))

foveal_radius_px = foveal_radius
foveal_radius_deg = foveal_radius * deg_per_pixel 

foveal_mask = np.zeros((peri_height,peri_width))

mid_point_w = np.round(peri_width/2)
mid_point_h = np.round(peri_height/2)

for i in range(0,peri_width):
	for j in range(0,peri_height):
		if np.sqrt(pow( i-mid_point_w,2) + pow(j-mid_point_h,2)) <= foveal_radius:
			foveal_mask[j][i] = 1

######################
# Create Filter bank #
######################

# To do.
# We have to add the gabor filter bank here.

visual_mask = 0

####################################
# Create list of Foveal_parameters #
####################################

foveal_param = foveal_param_init_all()

foveal_param.N_theta = N_theta;
foveal_param.N_e = N_e;
#foveal_param.wave_num = wave_num;
#foveal_param.orien = orien;
foveal_param.peri_height = peri_height;
foveal_param.peri_width = peri_width;

foveal_param.select_mask_stream = select_mask_stream;
foveal_param.select_mask_label = select_mask_label;
foveal_param.foveal_mask = foveal_mask;
foveal_param.foveal_radius_px = foveal_radius_px;
foveal_param.foveal_radius_deg = foveal_radius_deg;
#foveal_param.filter_bank = filter_bank;
#foveal_param.filter_bank_stream = filter_bank_stream;
foveal_param.peripheral_filters = peripheral_filters;

foveal_param.visual_mask = visual_mask;

foveal_param.deg_per_pixel = deg_per_pixel;

foveal_param.monitor.pixel_res_width = param.pixel_res_width;
foveal_param.monitor.pixel_res_height = param.pixel_res_height;
foveal_param.monitor.mon_width = param.mon_width;
foveal_param.monitor.mon_height = param.mon_height;
foveal_param.monitor.view_dist = param.view_dist;
