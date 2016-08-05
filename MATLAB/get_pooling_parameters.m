function [N_e N_theta] = get_pooling_parameters(scale,e0_in_deg,visual_field_radius_in_deg,deg_per_pixel)

w_theta = scale/2;
%N_theta = ceil(2*pi/w_theta); %Previous implementation
N_theta = floor(2*pi/w_theta);


e_0 = e0_in_deg;
visual_field_width = round(2*(visual_field_radius_in_deg./deg_per_pixel));
%e_r = visual_field_width*sqrt(2)/2*deg_per_pixel; %Use this option if you will not use the 
%create_regions_vector_function_smooth_FS instead of the
%create_regions_vector_function_smooth function.
e_r = visual_field_width/2*deg_per_pixel;


w_ecc = scale;
N_e = ceil((log(e_r) - log(e_0))/w_ecc); %Previous implementation
%N_e = floor((log(e_r) - log(e_0))/w_ecc);

