function regions = create_regions_vector_function_smooth_FS(e0_in_deg,e_max,visual_field_width,deg_per_pixel,N_theta,N_e)

%Comment for final function version.
%clc;close all;clear all;

%Parameters for testing:
%deg_per_pixel = 0.022;
%visual_field_radius_in_deg = 5; 
%e0_in_deg = .49;
%N_theta = 20;
%N_e = 4; 

%visual_field_radius_in_deg = e_max;

visual = 1;

%visual_field_width = round(2*(visual_field_radius_in_deg./deg_per_pixel));

center_r = round(visual_field_width/2);
center_c = center_r;

regions = zeros(visual_field_width,visual_field_width, N_theta, N_e);

visual_field_width_half = round((visual_field_width/2));
%visual_field_max = sqrt(2*visual_field_width_half^2);
visual_field_max = 2*visual_field_width_half;


x_vec = linspace(-1,1,visual_field_width);
y_vec = zeros(1,visual_field_width);

for i=1:visual_field_width
	if (x_vec(i)>=-0.75) && (x_vec(i)<-0.25)
		y_vec(i) = (cos(pi/2*((x_vec(i)+0.25)*2)))^2;
	elseif (x_vec(i)>=-0.25) && (x_vec(i)<0.25)
		y_vec(i) = 1;
	elseif (x_vec(i)>=0.25) && (x_vec(i)<0.75)
		y_vec(i) = 1-(cos(pi/2*((x_vec(i)-0.75)*2)))^2;
	else
		%Do nothing
	end
end

%Plot only one curve:
%N_theta = 20;
w_theta = 2*pi/N_theta;
t=1/2;
%t=0; %changed this since back to t=1/2 tiling was not smooth in center horizontal axis

h = zeros(1,visual_field_width);

arg_h = zeros(N_theta,visual_field_width+1+visual_field_width);
h_vec = zeros(N_theta,visual_field_width+1+visual_field_width);

%arg_h = zeros(N_theta,visual_field_width+1);%+1+visual_field_width);
%h_vec = zeros(N_theta,visual_field_width+1);%+1+visual_field_width);


%Check? Should be padded with NaN's?


%for j=0:(N_theta-1)
for j=0:N_theta
	for i=1:visual_field_width+1+visual_field_width
	%for i=1:visual_field_width
		%Going from 0 to 2*pi:
		arg_h(j+1,i) = (((i-1)*1.0/visual_field_width)*2*pi-((w_theta*j)+(w_theta*(1-t)/2)))/w_theta;
		%arg_h(j+1,i) = (((i)*1.0/visual_field_width)*2*pi-((w_theta*j)+(w_theta*(1-t)/2)))/w_theta-1;
	
		if arg_h(j+1,i)<-0.75
			h_vec(j+1,i) = 0;
		elseif (arg_h(j+1,i)>=-0.75) && (arg_h(j+1,i)<-0.25)
			h_vec(j+1,i) = (cos((pi/2)*((arg_h(j+1,i)+0.25)*2)))^2;
		elseif (arg_h(j+1,i)>=-0.25) && (arg_h(j+1,i)<0.25)
			h_vec(j+1,i) = 1;
		elseif (arg_h(j+1,i)>=0.25) && (arg_h(j+1,i)<0.75)
			h_vec(j+1,i) = 1-(cos((pi/2)*((arg_h(j+1,i)-0.75)*2)))^2;
		elseif arg_h(j+1,i)>0.75
			h_vec(j+1,i) = 0;
		else
			%Do nothing;
		end
	end
end

%Optional
%h_vec = h_vec(:,18:end);

if visual 
	figure();
	hold on;
	for i=1:N_theta
		plot(h_vec(i,1:visual_field_width),'LineWidth',2);
	end
	axis([0 visual_field_width 0 1]);
	title('Filter Value vs Polar angle');
	set(gca,'XTick',[0:visual_field_width/8:visual_field_width]);
end

%Now plot the Filter value of g(n) vs Retinal eccentricity

%e_0 = 1.5;
e_0 = e0_in_deg;
%e_r = N_theta;
%e_r = visual_field_width*sqrt(2)/2*deg_per_pixel;
e_r = visual_field_width/2*deg_per_pixel;

%20 degrees of visual angle:

N_ecc = N_e;
%N_ecc = 4;
w_ecc = (log(e_r) - log(e_0))/N_ecc;


arg_g = zeros(N_ecc,visual_field_width+1+visual_field_width);
g_vec = zeros(N_ecc,visual_field_width+1+visual_field_width);

%arg_g = zeros(N_ecc,visual_field_width+1);%+1+visual_field_width);
%g_vec = zeros(N_ecc,visual_field_width+1);%+1+visual_field_width);

%Check? should be padded with NaN's?


for j=0:(N_ecc-1)
%for j=0:N_ecc
	for i=1:visual_field_width+1+visual_field_width
	%for i=1:visual_field_width
		%Going from 0 to 2*pi:
		%arg_g(j+1,i) = (log((i-1)*e_r/(visual_field_width*sqrt(2)/2))-(log(e_0)+w_ecc*(j+1)))/w_ecc;
		%arg_g(j+1,i) = (log((i-1)*e_r/visual_field_width/sqrt(2))-(log(e_0)+w_ecc*(j)))/w_ecc;
		%arg_g(j+1,i) = (log((i-1)*e_r/(visual_field_width*sqrt(2)/2))-(log(e_0)+w_ecc*(j+1)))/w_ecc;
		%arg_g(j+1,i) = (log((i-1)*e_r/(visual_field_width*sqrt(2)))-(log(e_0)+w_ecc*(j+1)))/w_ecc;
		arg_g(j+1,i) = (log((i-1)*e_r/(visual_field_width))-(log(e_0)+w_ecc*(j+1)))/w_ecc;
	

		if arg_g(j+1,i)<-0.75
			g_vec(j+1,i) = 0;
		elseif (arg_g(j+1,i)>=-0.75) && (arg_g(j+1,i)<-0.25)
			g_vec(j+1,i) = (cos((pi/2)*((arg_g(j+1,i)+0.25)*2)))^2;
		elseif (arg_g(j+1,i)>=-0.25) && (arg_g(j+1,i)<0.25)
			g_vec(j+1,i) = 1;
		elseif (arg_g(j+1,i)>=0.25) && (arg_g(j+1,i)<0.75)
			g_vec(j+1,i) = 1-(cos((pi/2)*((arg_g(j+1,i)-0.75)*2)))^2;
		elseif arg_g(j+1,i)>=0.75
			g_vec(j+1,i) = 0;
		else
			%Do nothing;
		end
	end
end

if visual

	figure();
	hold on;
	for i=1:N_ecc
	%for i=1:N_ecc
		plot(g_vec(i,1:visual_field_width),'LineWidth',2);
	end
	axis([0 visual_field_width 0 1]);
	title('Filter Value vs Retinal Eccentricity');
	set(gca,'XTick',0:visual_field_width/8:visual_field_width);
end


h_vec = h_vec(:,1:end-1);
g_vec = g_vec(:,1:end-1);


%Now get the x,y coordinates from the polar coordinates

map = double(zeros(visual_field_width,visual_field_width));
map2 = double(zeros(visual_field_width,visual_field_width));


theta_temp_matrix = double(zeros(visual_field_width,visual_field_width));
ang_sign = +1;

map_hybrid = zeros(visual_field_width,visual_field_width,N_theta,N_ecc);
map_hybrid2 = zeros(visual_field_width,visual_field_width,N_theta,N_ecc);


ang_hybrid = zeros(visual_field_width,visual_field_width);
ecc_hybrid = zeros(visual_field_width,visual_field_width);

for i=1:visual_field_width
	for j=1:visual_field_width

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Optimized implementation %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		%Get distance from center
		dist = sqrt((visual_field_width_half-i)^2+(visual_field_width_half-j)^2);

		%get angle from center
		if i~=visual_field_width_half && j~=visual_field_width_half

			ang = atan2(visual_field_width_half-i,j-visual_field_width_half);
			true_ang = ang;

			if i<visual_field_width_half && j<visual_field_width_half
				true_ang = pi + ang;
			end

			if i>=visual_field_width_half && j>=visual_field_width_half
				true_ang = pi + ang;
			end

			if i<=visual_field_width_half && j>visual_field_width_half
				true_ang = pi + ang;
			end
		
			if i>visual_field_width_half && j<=visual_field_width_half
				true_ang = pi + ang;
			end

		else
			ang = 0;
			true_ang = ang;
		end

		if i<=visual_field_width_half && j==visual_field_width_half
			true_ang = atan2(visual_field_width_half-i,1) + pi;
		end

		if j>visual_field_width_half && i==visual_field_width_half
			true_ang = atan2(visual_field_width_half-i,1) + pi;
		end

		if j==visual_field_width_half && i>visual_field_width_half
			true_ang = atan2(visual_field_width_half-i,1) + pi;
		end


		%find closest angle match:
		ang_match = round(true_ang/(2*pi)*visual_field_width);

		%find closest eccentricity match:
		dist_match = round(dist/(visual_field_width/2)*visual_field_width);


		if ang_match<=0
			ang_match = 1;
		end

		if dist_match<=0
			dist_match = 1;
		end

		%Get Hybrid Computations
		ang_hybrid(i,j) = ceil(true_ang/(2*pi)*N_theta);

		%if ang_hybrid(i,j) == 0
		if ang_hybrid(i,j) <= 0
			ang_hybrid(i,j) = 1;
		end

		ecc_hybrid(i,j) = ceil(dist_match/visual_field_max/2*N_ecc);

		%Added in Version 4
		if ecc_hybrid(i,j)>N_ecc
			ecc_hybrid(i,j) = N_ecc;
		end

		%Also find theta value and eccentricity
		theta_temp = ceil(ang_match/(visual_field_width/N_theta));
		true_ang_matrix(i,j) = true_ang;

		ang_theta_matrix(i,j) = floor(true_ang*N_theta/(2*pi))+1;
		
		temp_ecc = find(max(g_vec(:,dist_match))==g_vec(:,dist_match));

		if length(temp_ecc)>1

			ecc_matrix(i,j) = 1;
		else
			ecc_matrix(i,j) = temp_ecc;
		end

		if ecc_matrix(i,j)==(N_theta+1)
			ecc_matrix(i,j) = N_theta;
		end

		temp_theta = find(max(h_vec(:,ang_match))==h_vec(:,ang_match));

		if length(temp_theta)>1
			theta_matrix(i,j) = 1;
		else
			theta_matrix(i,j) = temp_theta;
		end

		if theta_matrix(i,j) == (N_theta+1)
			theta_matrix(i,j) = 1;
		end

		h_buffer_indx = find(h_vec(:,ang_match)>0);
		g_buffer_indx = find(g_vec(:,dist_match)>0);

		% Get Shadow tones of every region. This is useful for Metamer Construction
		% Since every pooling region is locally constructed and then blended
			
		for z1 = 1:length(h_buffer_indx)
			for z2 = 1:length(g_buffer_indx)
				map_hybrid2(i,j,h_buffer_indx(z1),g_buffer_indx(z2)) = h_vec(h_buffer_indx(z1),ang_match).*g_vec(g_buffer_indx(z2),dist_match);

			end
		end

		hybrid_buffer = map_hybrid2(i,j,:,:);
		map2(i,j) = max(hybrid_buffer(:));
	end
end

if visual 
	figure();
	imshow(map2);

	figure();
	imagesc(map2), colorbar;
end


%Final Output:
regions = map_hybrid2;

