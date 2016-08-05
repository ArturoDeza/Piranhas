-- Create_Pianhas.lua script:

-- Load the piranhas module
require 'piranhas'

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

create_regions_vector_smooth(e0_in_deg,e_max,visual_field_width,deg_per_pixel,N_theta,N_e)
