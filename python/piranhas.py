import random
import numpy
import math
import numpy as np
import matplotlib.pyplot as plt

# Define monitor class:
class monitor:
	# check pixel resolution of the computer where
	# phsychophysics will be used.
	pixel_res_width = None
	pixel_res_height = None

	# specify these values in cm of the visible
	# screen area of the monitor
	mon_width = None
	mon_height = None 
	view_dist = None

	def __init__(self,pixel_res_width,pixel_res_height,mon_width,mon_height,view_dist):
		self.pixel_res_width = pixel_res_width
		self.pixel_res_height = pixel_res_height
		self.mon_width = mon_width
		self.mon_height = mon_height
		self.view_dist = view_dist

# Define param class:
class param:
	# Set visual parameters to zero by default
	visual = 0
	visual_mask = 0

	# Default settings
	e0_in_deg = 0.25

	# You should initialize the monitor first
	# to avoid any errors
	cm_per_pixel = None
	deg_per_pixel = None

	def __init__(self,monitor_input,visual_field_radius_in_deg,fovea,scale):
		# Initialize monitor parameters:
		self.monitor = monitor_input

		#Initialize Field of View parameters:
		self.visual_field_radius_in_deg = visual_field_radius_in_deg
		self.fovea = fovea
		self.scale = scale
		self.cm_per_pixel = monitor_input.mon_width/monitor_input.pixel_res_width
		self.deg_per_pixel = 2*math.degrees(math.atan(self.cm_per_pixel/2/monitor_input.view_dist))

	def sample_V1(self):
		# Field of View paramaters:
		self.visual_field_radius_in_deg = 10
		self.fovea = 2.0
		self.scale = 0.25

	def sample_V2(self):
		# Field of View paramaters:
		self.visual_field_radius_in_deg = 10
		self.fovea = 2.0
		self.scale = 0.25

	def set_fovea_deg(self,foveal_input):
		self.fovea = foveal_input

	def set_fovea_px(self,foveal_input):
		self.fovea = foveal_input/deg_per_pixel

	def set_scale(self,scale_input):
		# Try to make these values to be between 
		# 0.25 to 0.5 as a good rule of thumb.
		self.scale = scale_input

	def get_scale(self):
		return self.scale

	# Still haven't tested how good this function works
	# After pooling regions inside fovea are removed (which 
	# should be none).
	def nullify_e0_deg(self):
		e0_in_deg = self.fovea

	def set_e0_deg(self,e0_in_deg_input):
		e0_in_deg = e0_in_deg_input

	# Turn visuals on/off:
	def visual_on(self):
		visual = 1

	def visual_off(self):
		visual = 0
	
	#Turn visual_masks on/off:
	def visual_mask_on(self):
		visual_mask = 1

	def visual_mask_off(self):
		visual_mask = 0

# Define Foveal Param class:
class foveal_param:
	N_theta = None
	N_e = None
	peri_height = None
	peri_width = None

	select_mask_stream = None #cell
	select_mask_label = None #cell
	foveal_mask = None
	foveal_radius_px = None
	foveal_radius_deg = None

	peripheral_filters = None #cell
	visual_mask = 0
	deg_per_pixel = None


	# All of these parameters will be passed after the initial params are defined and computed
	def __init__(self,monitor_input,N_theta,N_e,peri_height,peri_width,select_mask_stream,select_mask_label,foveal_mask,foveal_radius_px,foveal_radius_deg):
		self.N_theta = N_theta
		self.N_e = N_e
		self.peri_height = peri_height
		self.peri_width = peri_width
		self.select_mask_stream = select_mask_stream
		self.select_mask_label = select_mask_label
		self.foveal_mask = foveal_mask
		self.foveal_radius_px = foveal_radius_px
		self.foveal_radius_deg = foveal_radius_deg
		self.monitor_settings = monitor_input

		self.cm_per_pixel = monitor_input.mon_width/monitor_input.pixel_res_width
		self.deg_per_pixel = 2*math.degrees(math.atan(self.cm_per_pixel/2/monitor_input.view_dist))
	
###########################################
# Compute Eccentricity + Theta parameters #
###########################################

def get_pooling_parameters(scale,e0_in_deg,visual_field_radius_in_deg,deg_per_pixel):
  
  w_theta = scale/2
  N_theta = math.floor(2*np.pi/w_theta)
  
  e_0 = e0_in_deg
  visual_field_width = round(2*(visual_field_radius_in_deg/deg_per_pixel))
  e_r = visual_field_width/2*deg_per_pixel
  
  w_ecc = scale
  N_e = math.ceil((math.log(e_r) - math.log(e_0))/w_ecc) # Holding this for the previous implementation
  
  return (N_e, N_theta)

###########################
# Create Regions function #
###########################

def create_regions_vector_function_smooth(e0_in_deg,e_max,visual_field_width,deg_per_pixel,N_theta,N_e):

	# visual flag on?

	visual = 1

	center_r = round(visual_field_width/2)
	center_c = center_r

	regions = np.zeros((visual_field_width,visual_field_width,N_theta,N_e))

	visual_field_width_half = round(visual_field_width*1.0/2) # before this was not multiplied by 1.0
	visual_field_max = np.sqrt(2*pow(visual_field_width_half,2))

	x_vec = np.linspace(-1,1,visual_field_width)
	y_vec = np.zeros(visual_field_width)


	for i in range(1,int(round(visual_field_width))):
		if (x_vec[i] >=-0.75) and (x_vec[i] < -0.25):
			y_vec[i] = pow(math.cos((np.pi/2)*((x_vec[i]+0.25)*2)),2)
		elif (x_vec[i] >= -0.25) and (x_vec[i] < 0.25):
			y_vec[i] = 1.0
		elif (x_vec[i] >= 0.25) and (x_vec[i] < 0.75):
			y_vec[i] = -pow(math.cos((np.pi/2)*((x_vec[i]+0.75)*2)),2)

	# Initialize hyperparameters for h function

	w_theta = 2*np.pi/N_theta
	t = 0

	# Creating the h function.

	h = np.zeros(visual_field_width)

	arg_h = np.zeros((N_theta,visual_field_width + 1 + visual_field_width))
	h_vec = np.zeros((N_theta,visual_field_width + 1 + visual_field_width))

	for j in range(0,int(round(N_theta-1))+1):
		for i in range(1,int(round(visual_field_width + 1 + visual_field_width))):
		
			arg_h[j,i] = (((i-1)*1.0/visual_field_width)*2*np.pi - ((w_theta*j)+(w_theta*(1-t)/2)))/w_theta
		
			if arg_h[j,i]<-0.75:
				h_vec[j,i] = 0
			elif arg_h[j,i]>=-0.75 and arg_h[j,i]<-0.25:
				h_vec[j,i] = pow(math.cos((np.pi/2)*((arg_h[j,i]+0.25)*2)),2)
			elif arg_h[j,i]>=-0.25 and arg_h[j,i]<0.25:
				h_vec[j,i] = 1
			elif arg_h[j,i]>=0.25 and arg_h[j,i]<0.75:
				h_vec[j,i] = 1 - pow(math.cos((np.pi/2)*((arg_h[j,i]-0.75)*2)),2)
			elif arg_h[j,i]>0.75:
				h_vec[j,i] = 0

	# start doing visualization of h_vec piranha functions:

	if visual:
		for i in range(0,int(round(N_theta-1))+1):
			plt.plot(h_vec[i,:])
		plt.show()

	# Initialize hyperparameters for g function

	e_0 = e0_in_deg
	e_r = visual_field_width*math.sqrt(2)/2*deg_per_pixel

	N_ecc = N_e
	w_ecc = (math.log(e_r)-math.log(e_0))/N_ecc

	arg_g = np.zeros((N_ecc,visual_field_width+1+visual_field_width))
	g_vec = np.zeros((N_ecc,visual_field_width+1+visual_field_width))

	for j in range(0,int(round(N_ecc-1))+1):
		for i in range(1, int(round(visual_field_width + 1 + visual_field_width))):
		
			arg_g[j,i] = (math.log((i-1+0.00001)*e_r/(visual_field_width*math.sqrt(2))) - (math.log(e_0)+w_ecc*(j)))/w_ecc
		
			if arg_g[j,i] < -0.75:
				g_vec[j,i] = 0
			elif arg_g[j,i] >= -0.75 and arg_g[j,i]<-0.25:
				g_vec[j,i] = pow(math.cos((np.pi/2)*((arg_g[j,i]+0.25)*2)),2)
			elif arg_g[j,i] >= -0.25 and arg_g[j,i] < 0.25:
				g_vec[j,i] = 1
			elif arg_g[j,i] >= 0.25 and arg_g[j,i] < 0.75:
				g_vec[j,i] = 1 - pow(math.cos((np.pi/2)*((arg_g[j,i]-0.75)*2)),2)
			elif arg_g[j,i] >= 0.75:
				g_vec[j,i] = 0

	# Add visualization map:

	if visual:
		for i in range(0,int(round(N_ecc-1))+1):
			plt.plot(g_vec[i,:])
		plt.show()

	##########################################################
	# Now get the x,y coordinates from the polar coordinates #
	##########################################################

	map = np.zeros((visual_field_width,visual_field_width))
	map2 = np.zeros((visual_field_width,visual_field_width))

	theta_temp_matrix = np.zeros((visual_field_width,visual_field_width))

	ang_sign = 1

	#Create empty maps:

	map_hybrid = np.zeros((((visual_field_width,visual_field_width,N_theta,N_ecc))))
	map_hybrid2 = np.zeros((((visual_field_width,visual_field_width,N_theta,N_ecc))))

	ang_hybrid = np.zeros((visual_field_width,visual_field_width))
	ecc_hybrid = np.zeros((visual_field_width,visual_field_width))

	true_ang_matrix = np.zeros((visual_field_width,visual_field_width))
	ang_theta_matrix = np.zeros((visual_field_width,visual_field_width))

	ecc_matrix = np.zeros((visual_field_width,visual_field_width))
	theta_matrix = np.zeros((visual_field_width,visual_field_width))

	for i in range(0,int(round(visual_field_width))):
		for j in range(0,int(round(visual_field_width))):
		
			#############################
			# Optimized Implementation: #
			#############################
		
			# Get Distance from center:
			dist = math.sqrt(pow(visual_field_width_half-i,2)+pow(visual_field_width_half-j,2))
		
			# Get angle from center:
			if i!=visual_field_width_half and j!=visual_field_width_half:
			
				ang = math.atan2(visual_field_width_half-i,j-visual_field_width_half)
				true_ang = ang
			
				if i<visual_field_width_half and j < visual_field_width_half:
					true_ang = np.pi + ang
			
				if i>=visual_field_width_half and j>= visual_field_width_half:
					true_ang = np.pi + ang
						
				if i<=visual_field_width_half and j >visual_field_width_half:
					true_ang = np.pi + ang
			
				if i>visual_field_width_half and j<= visual_field_width_half:
					true_ang = np.pi + ang
		
			else:
				ang = 0
				true_ang = ang
		
			if i<=visual_field_width_half and j == visual_field_width_half:
				true_ang = math.atan2(visual_field_width_half-i,1) + np.pi
		
			if j>visual_field_width_half and i == visual_field_width_half:
				true_ang = math.atan2(visual_field_width_half-i,1) + np.pi
		
			if j == visual_field_width_half and i > visual_field_width_half:
				true_ang =math.atan2(visual_field_width_half-i,1) + np.pi
		
			# Find closest angle match:
			ang_match = round(true_ang/(2*np.pi)*visual_field_width)
		
			#Find closest eccentricty match:
			dist_match = round(dist/(visual_field_width/2)*visual_field_width)
		
			if ang_match<=0:
				ang_match = 1
		
			if dist_match<=0:
				dist_match = 1
		
			# Get Hybrid Computations:
			ang_hybrid[i][j] = math.ceil(true_ang/(2*np.pi)*N_theta)
		
			if ang_hybrid[i][j] <= 0:
				ang_hybrid[i][j] = 1
		
			ecc_hybrid[i][j] = math.ceil(dist_match/visual_field_max/2*N_ecc)
		
			if ecc_hybrid[i][j]>N_ecc:
				ecc_hybrid[i][j] = N_ecc
		
			# Also find theta value and eccentricity:
			theta_temp = math.ceil(ang_match/(visual_field_width/N_theta))
			true_ang_matrix[i][j] = true_ang
		
			ang_theta_matrix[i][j] = math.floor(true_ang*N_theta/(2*np.pi))+1
		
			temp_ecc = np.where(max(g_vec[:,dist_match])==g_vec[:,dist_match])
		
			if len(temp_ecc[0]) > 1:
				ecc_matrix[i][j] = 0
			else:
				ecc_matrix[i][j] = temp_ecc[0]
		
			if ecc_matrix[i][j] == (N_ecc+1):
				ecc_matrix[i][j] = N_ecc
		
			temp_theta = np.where(max(h_vec[:,ang_match])==h_vec[:,ang_match])
		
			if len(temp_theta[0]) > 1:
				theta_matrix[i][j] = 0
			else:
				theta_matrix[i][j] = temp_theta[0]
		
			if theta_matrix[i][j] == (N_theta + 1):
				theta_matrix[i][j] = 0
		
			h_buffer_indx = np.where(h_vec[:,ang_match]>0)
			g_buffer_indx = np.where(g_vec[:,dist_match]>0)
		
			# Get Smooth shadow tones for every regions::
			for z1 in range(0,len(h_buffer_indx[0])):
				for z2 in range(0,len(g_buffer_indx[0])):
					map_hybrid2[i][j][h_buffer_indx[0][z1]][g_buffer_indx[0][z2]] = h_vec[h_buffer_indx[0][z1],ang_match] * g_vec[g_buffer_indx[0][z2],dist_match]
		
			hybrid_buffer = map_hybrid2[i][j][:][:]
			map2[i][j] = max(hybrid_buffer.flatten())

	if visual:
		plt.imshow(ecc_matrix)
		plt.show()

	if visual:
		plt.imshow(theta_matrix)
		plt.show()

	if visual:
		plt.imshow(map2)
		plt.show()
	
	regions = map_hybrid2

	return regions

##################################################
# Generate pooling regions data structure for    #
# multiple operations. Excludes regions          #
# inside the foveal radius.                      #
##################################################

def generate_pooling_regios_vector_smooth(deg_per_pixel,N_e,N_theta,visual_field_radius_in_deg,fovea,e0_in_deg,visual):

	e_max = visual_field_radius_in_deg
	visual_field_width = round(2*(visual_field_radius_in_deg/deg_per_pixel))

	center_r = round(visual_field_width/2)
	center_c = center_r

	regions = create_regions_vector_function_smooth(e0_in_deg,e_max,visual_field_width,deg_per_pixel,N_theta,N_e)

	regions = regions.transpose([2,3,0,1])
	#fovea = 1.0 #uncomment to pilot studies

	centers = np.zeros(((2,N_theta,N_e)))
	areas = np.zeros((N_theta,N_e))

	mask_matrix = np.zeros((N_theta,N_e))

	mask = np.zeros((visual_field_width,visual_field_width))

	#Double check this:

	for nt in range(0,int(round(N_theta))):
		for ne in range(0,int(round(N_e))):
			#Verify if this squeeze is equivalent to MATLAB's squeeze operation:
			mask  = regions[nt][ne][:][:].squeeze()
			idx_rc = np.where(mask)
			r = idx_rc[0]
			c = idx_rc[1]
	
			centers[0,nt,ne] = np.mean(r)
			centers[1,nt,ne] = np.mean(c)
	
			areas[nt,ne] = len(r)

	class Filters:
		def __init__(self,N_theta,N_e,visual_field_width):
			self.regions = np.array((((N_theta,N_e,visual_field_width,visual_field_width))))
			self.centers = np.array(((2,N_theta,N_e)))
			self.areas = np.array((N_theta,N_e))
			self.offsets_r = [[None for j in range(0,int(round(N_e)))] for i in range(0,int(round(N_theta)))] #empty list
			self.offsets_c = [[None for j in range(0,int(round(N_e)))] for i in range(0,int(round(N_theta)))] #empty list
			self.weights = [[None for j in range(0,int(round(N_e)))] for i in range(0,int(round(N_theta)))] #empty list
			self.uniq_pix = [[None for j in range(0,int(round(N_e)))] for i in range(0,int(round(N_theta)))] #empty list

	filters = Filters(N_theta,N_e,visual_field_width)

	filters.regions = regions
	filters.centers = centers
	filters.areas = areas

	# We want to preserve smooth curves:
	blindspot_threshold = 0.0

	# Offset coordinates for pooling regions:
	center_r = round(filters.regions.shape[2]/2)
	center_c = round(filters.regions.shape[3]/2)

	#Initialize filters_offsets, filter weights, and filters uniq_pix

	for nt in range(0,filters.regions.shape[0]):
		for ne in range(0,filters.regions.shape[1]):
			rc_buff = np.where(np.squeeze(filters.regions[nt,ne,:,:])>0)
	
			rs = rc_buff[0]
			cs = rc_buff[1]
	
			filters.offsets_r[nt][ne] = rs - center_r
			filters.offsets_c[nt][ne] = cs - center_c
	
			idx1 = nt*np.ones(rs.shape[0])
			idx2 = ne*np.ones(rs.shape[0])
			idx3 = rs
			idx4 = cs
	
			idx1 = idx1.astype(int)
			idx2 = idx2.astype(int)
			idx3 = idx3.astype(int)
			idx4 = idx4.astype(int)
	
			weights_buff = filters.regions[idx1,idx2,idx3,idx4]
			filters.weights[nt][ne] = weights_buff.flatten()
	
			rc_buff2 = np.where(np.squeeze(filters.regions[nt][ne][:][:]))
			filters.uniq_pix[nt][ne] = rc_buff2

	# Time to do some painting and paint the unique pixels:

	# Define empty colored map
	foveal_map_colored = np.zeros((int(visual_field_width),int(visual_field_width),3),dtype = np.uint8)

	import matplotlib.pyplot as plt
	import matplotlib.cm as cmx
	import matplotlib.colors as colors

	def get_cmap(N):
		color_norm = colors.Normalize(vmin=0,vmax=N-1)
		scalar_map = cmx.ScalarMappable(norm=color_norm,cmap='hsv')
		def map_index_to_rgb_color(index):
			return scalar_map.to_rgba(index)
		return map_index_to_rgb_color


	cmap = get_cmap(filters.regions.shape[0] * filters.regions.shape[1])

	Nt_range = filters.regions.shape[0]
	Ne_range = filters.regions.shape[1]

	#for nt in range(0,Nt_range):
	#	for ne in range(0,Ne_range):
	#		for z in range(0,filters.uniq_pix[nt][ne][0].shape[0]):
	#			foveal_map_colored[filters.uniq_pix[nt][ne][0][z],filters.uniq_pix[nt][ne][1][z],0] = 255*cmap((ne-1)*filters.regions.shape[1]+nt)[0]
	#			foveal_map_colored[filters.uniq_pix[nt][ne][0][z],filters.uniq_pix[nt][ne][1][z],1] = 255*cmap((ne-1)*filters.regions.shape[1]+nt)[1]
	#			foveal_map_colored[filters.uniq_pix[nt][ne][0][z],filters.uniq_pix[nt][ne][1][z],2] = 255*cmap((ne-1)*filters.regions.shape[1]+nt)[2]

	#from PIL import Image
	#
	#if visual == 1:
	#	img = Image.fromarray(foveal_map_colored,'RGB')
	#	img.show()

	#Some Lines missing from MATLAB code.
	#Will add them when relevant.

	peripheral_filters = filters
	fovea_radius = fovea

	# Decide up to which cell to discard (along the radial axis)
	for n_e in range(0,Ne_range):
		offsets_r = abs(peripheral_filters.offsets_r[0][n_e])
		offsets_c = abs(peripheral_filters.offsets_c[0][n_e])
		# See how much of this cell is within the fovea:
		within = (offsets_r<=fovea_radius/deg_per_pixel) & (offsets_c<=fovea_radius/deg_per_pixel)
	
		if sum(within)/len(within)<0.5:
			n_e_limit = n_e -1 #subtracing 1 here since we have 0 indexing.
			break

	#Now use the value from n_e to discard pooling regions (not including): n_e
	#Essentially changing the dimension of peripheral_filters.
	#Might have to create a new data structure here.


	new_extension = Ne_range - n_e + 1

	peripheral_filters.regions = peripheral_filters.regions[:,n_e_limit:,:,:]
	peripheral_filters.centers = peripheral_filters.centers[:,:,n_e_limit:]
	peripheral_filters.areas = peripheral_filters.areas[:,n_e_limit:]

	# Have to do this since python 2.7 doesn't support MultiDimensional List Slicing
	# If someone has a workaround go for it.

	offsets_r_buff = [[None for j in range(0,int(round(new_extension)))] for i in range(0,int(round(N_theta)))] 
	offsets_c_buff = [[None for j in range(0,int(round(new_extension)))] for i in range(0,int(round(N_theta)))] 
	weights_buff = [[None for j in range(0,int(round(new_extension)))] for i in range(0,int(round(N_theta)))] 
	uniq_pix_buff = [[None for j in range(0,int(round(new_extension)))] for i in range(0,int(round(N_theta)))] 

	offsets_r_buff = [[peripheral_filters.offsets_r[i][j] for j in range(n_e_limit,int(N_e))] for i in range(0,int(N_theta))]
	offsets_c_buff = [[peripheral_filters.offsets_c[i][j] for j in range(n_e_limit,int(N_e))] for i in range(0,int(N_theta))]
	weights_buff = [[peripheral_filters.weights[i][j] for j in range(n_e_limit,int(N_e))] for i in range(0,int(N_theta))]
	uniq_pix_buff = [[peripheral_filters.uniq_pix[i][j] for j in range(n_e_limit,int(N_e))] for i in range(0,int(N_theta))]


	peripheral_filters.offsets_r = offsets_r_buff
	peripheral_filters.offsets_c = offsets_c_buff
	peripheral_filters.weights = weights_buff
	peripheral_filters.uniq_pix = uniq_pix_buff

	#discard low weight pixels
	weight_threshold = 0.0

	set_zero_indx = np.where(peripheral_filters.regions<=weight_threshold)
	peripheral_filters.regions[set_zero_indx[0],set_zero_indx[1],set_zero_indx[2],set_zero_indx[3]] = 0.0 

	offset_new_range = np.shape(peripheral_filters.offsets_r)

	# All pixels should be valid:

	# Make sure this last block is correct

	for i in range(0,offset_new_range[0]):
		for j in range(0,offset_new_range[1]):
	
			valids = peripheral_filters.weights[i][j]>weight_threshold
			valids_indx = np.where(valids)
	
			#Reassignment begins here:		
			peripheral_filters.offsets_r[i][j] = peripheral_filters.offsets_r[i][j][valids_indx]
			peripheral_filters.offsets_c[i][j] = peripheral_filters.offsets_c[i][j][valids_indx]
			peripheral_filters.weights[i][j] = peripheral_filters.weights[i][j][valids_indx]
			peripheral_filters.areas[i,j] = sum(valids)		

	# Have the extra visualization of the suppressed regions inside fovea.
	#
	#if visual:
	#	visual_field_width = np.shape(peripheral_filters.regions,3)

	num_elements = offset_new_range[0]*offset_new_range[1]

	W = np.reshape(peripheral_filters.regions,(num_elements,visual_field_width,visual_field_width),order="F")
	W2 = np.max(W,0)

	plt.imshow(W2)
	plt.show()


	return peripheral_filters
