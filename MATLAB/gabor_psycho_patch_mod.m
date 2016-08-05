%function [g_mod_patch g_real_patch g_imag_patch] = gabor_psycho_patch(mon_width,mon_height,pixel_res_width,pixel_res_height,view_dist,cycles_per_deg,theta_c,psi_c,bandwidth_deg,gamma_c)

function [g_mod_patch g_real_patch g_imag_patch] = gabor_psycho_patch_mod(pixel_res_width,pixel_res_height,mon_width,mon_height,view_dist,cycles_per_deg_wave,bandwidth,gamma_c,theta_c,psi_c)
% Function written by Arturo Deza @ VIU-LAB, UCSB
% This function is targeted to creating Gabor patched for psychophysics experiments

% The following units should be specified in (cm), 
% but could also be specified in another units like inches. 
% The note to keep in mind is that all length units should be consistent.

% mon_width = monitor_width that you will be using to should psychophysical stimuli; (cm)
% mon_height = monitor_height that you will be using to should psychophysical stimuli; (cm)
% view_dist = viewing_distance between the monitor and the participants eyes -- that you will be using to should psychophysical stimuli; (cm)
% cycles_per_degree: cycles per degree of the sine modulation Gabor patch


% lambda_c = wavelength (in pixels);
% theta_c = Orientation with respect to vertical;
% psi_c = phase_offset (in radians);
% sigma_c = standard deviation of the gaussian envelope (in pixels);
% gamma_c = spatial aspect ratio (usually want this value to be 1 or sometimes 2)

% create deg_per_pixel conversion
% You do not need the width for the deg_per_pixel modulation,
% unless you want to modulate resolution by width
deg_per_pixel = 2*atand(mon_height/2/view_dist)/pixel_res_height;
deg_per_pixel_w = 2*atand(mon_width/2/view_dist)/pixel_res_width;

%if deg_per_pixel resolution has a greater difference than tolerance, this should be a problem!

deg_pix_tolerance = 0.05;
if deg_per_pixel-deg_per_pixel_w>=deg_pix_tolerance
	warning('Check pixel and monitor parameters!');
end

lambda_c = 1/(cycles_per_deg_wave*deg_per_pixel);

% Old style bandwidth in degrees
%sigma_c = 1/(2*deg_per_pixel)*bandwidth_deg;

% USing bandwidth in octaves
%sigma_c = lambda_c * sqrt(log(2)/2)*(2^band_octave_envelope+1)/(2^band_octave_envelope-1)/pi;

sigma_c = lambda_c*sqrt(log(2)/2)*(2^bandwidth+1)/(2^bandwidth-1)/pi;


%disp(lambda_c);
%disp(sigma_c);

%Since the value of a gaussian after 5 sigma's is extremly low:
w_patch = round(sigma_c*5);
h_patch = round(sigma_c*5);

%g_function
g_real_patch = zeros(w_patch,h_patch);
g_imag_patch = zeros(w_patch,h_patch);
g_mod_patch = zeros(w_patch,h_patch);

% Comment this out for a bit to see if vectorized solution works

for i=1:w_patch %this will be the x proxy value	
	for j=1:h_patch %this will be the y proxy value

		x_p = ((i-w_patch/2)*cos(theta_c) + (j-h_patch/2)*sin(theta_c));
		y_p = (-(i-w_patch/2)*sin(theta_c) + (j-h_patch/2)*cos(theta_c));

		g_real_patch(i,j) = exp(-(x_p^2+gamma_c^2*y_p^2)/(2*sigma_c^2))*cos(2*pi*x_p/lambda_c + psi_c);	
		g_imag_patch(i,j) = exp(-(x_p^2+gamma_c^2*y_p^2)/(2*sigma_c^2))*sin(2*pi*x_p/lambda_c + psi_c);

	end
end

%Vectorized form
%[x_p2 y_p2] = meshgrid(-round(w_patch/2):1:round(w_patch/2),-round(h_patch/2):1:round(h_patch/2));
%
%x_p = x_p2*cos(theta_c) + y_p2*sin(theta_c);
%y_p = -x_p2*sin(theta_c) + y_p2*cos(theta_c);
%
%g_real_patch = exp(-(x_p.^2+gamma_c^2*y_p.^2)/(2*sigma_c^2))*cos(2*pi*x_p/lambda_c + psi_c);
%g_imag_patch = exp(-(x_p.^2+gamma_c^2*y_p.^2)/(2*sigma_c^2))*sin(2*pi*x_p/lambda_c + psi_c);



g_real_patch = g_real_patch - mean(g_real_patch(:));%*ones(w_patch,h_patch);
g_imag_patch = g_imag_patch - mean(g_imag_patch(:));%*ones(w_patch,h_patch);

g_mod_patch = sqrt(g_real_patch.^2 + g_imag_patch.^2);
