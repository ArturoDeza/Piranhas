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
