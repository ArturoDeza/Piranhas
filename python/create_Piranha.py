# Python Script to initialize a Peripheral Architecture:
# Developed by Arturo Deza as part of the Piranhas Toolkit
# Questions & bugs: deza @ dyns.ucsb.edu

from piranhas import *

# Run Intialization Parameters:
param = param_init_all()

scale = param.scale
fovea = param.fovea
e0_in_deg = param.e0_in_deg
visual_field_radius_in_deg = param.visual_field_radius_in_deg
deg_per_pixel = param.deg_per_pixel


# Get Peripheral Architecture Parameters

N_param = get_pooling_parameters(scale,e0_in_deg,visual_field_radius_in_deg,deg_per_pixel)

N_e = N_param[0]
N_theta = N_param[1]

# Create Regions:

# Commenting this out since it went on the previous version
#
#e_max = visual_field_radius_in_deg
#visual_field_width = round(2*(visual_field_radius_in_deg/deg_per_pixel))
#
#center_r = round(visual_field_width/2)
#center_c = center_r
#
#regions = create_regions_vector_function_smooth(e0_in_deg,e_max,visual_field_width,deg_per_pixel,N_theta,N_e)

##############################
# Indexed Peripheral Regions #
##############################

visual = 0

peripheral_regions = generate_pooling_regios_vector_smooth(deg_per_pixel,N_e,N_theta,visual_field_radius_in_deg,fovea,e0_in_deg,visual)

# The peripheral_regions data structure allows you to do anything you want with the pooling regions.



