% This script will create a Gabor Filter Bank
%clc;close all;clear all;

%function filter_bank = create_gabor_bank_mod(orien,env_zero,oct_num_env,wave_zero,oct_num_wave,param,vis)
function filter_bank = create_gabor_bank_mod(orien,freq_zero,oct_freq,bandwidth,param,vis)
%vis = 1;
%oct_num = 4
%orien = 4;

%param.pixel_res_width = 1240;
%param.pixel_res_height = 1024;
%param.mon_width = 33;
%param.mon_height = 33;
%param.view_dist = 76;
%param.cycles_per_deg = 2;
%
%param.gamma_c = 1;
%
%param.psi_c =0;

oct_freq_vec = [freq_zero];

%number of octaves:
for i=0:oct_freq-2
	%oct_freq_vec = [oct_freq_vec; oct_freq_vec(i+1)*2^(i+1)];
	oct_freq_vec = [oct_freq_vec; freq_zero*2^(i+1)];
end


orien_vec = linspace(0,pi,orien+1);
orien_vec = orien_vec(1:end-1);

pixel_res_width = param.pixel_res_width;
pixel_res_height = param.pixel_res_height;
mon_width = param.mon_width;
mon_height = param.mon_height;
view_dist = param.view_dist;
%cycles_per_deg = param.cycles_per_deg;
%band_octave = param.band_octave;
gamma_c = param.gamma_c;
%theta_c = param.theta_c;
psi_c = param.psi_c;



for i=1:length(oct_freq_vec)
		for j=1:length(orien_vec)
			[g_mod g_real g_imag] = gabor_psycho_patch_mod(pixel_res_width,pixel_res_height,mon_width,mon_height,view_dist,oct_freq_vec(i),bandwidth,gamma_c,orien_vec(j),psi_c);
			filter_bank{i}{j}{1} = g_real;
			filter_bank{i}{j}{2} = g_imag; 
		end
end



if vis == 1
	for i=1:length(oct_freq_vec)
		for j=1:length(orien_vec)
			[size_bank_w size_bank_h] = size(filter_bank{i}{j}{1});

			bank_matrix{i}(1:size_bank_h,1+(j-1)*size_bank_w:j*size_bank_w) = filter_bank{i}{j}{1};
			bank_matrix{i}(size_bank_h+1:2*size_bank_h,1+(j-1)*size_bank_w:j*size_bank_w) = filter_bank{i}{j}{2};
		end

		scale_str = num2str(i);

		figure();imshow(mat2gray(bank_matrix{i}));
	end
end
