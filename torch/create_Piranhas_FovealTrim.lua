-- Create_Pianhas.lua script:

-- Load the piranhas module
require 'piranhas'
require 'torch'
require 'gnuplot'

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

regions = create_regions_vector_smooth(e0_in_deg,e_max,visual_field_width,deg_per_pixel,N_theta,N_e)

regions = regions:permute(3,4,1,2)
-- fovea = 1.0 -- Un comment for debugging purposes

for ne=1,N_e do
	regions[{1,ne,{},{}}] = regions[{1,ne,{},{}}] + regions[{N_theta+1,ne,{},{}}]
end	

-- Discard the supplementary region:
regions = regions[{{1,N_theta},{},{},{}}]

centers = torch.zeros(2,N_theta,N_e)
areas = torch.zeros(N_theta,N_e)

mask_matrix = torch.zeros(N_theta,N_e)

mask = torch.zeros(visual_field_width,visual_field_width)

-- Double Check this: [Adding scropes]
do
	local nt=1
	local ne=1
	for nt=1,N_theta do
		for ne=1,N_e do
		
			mask = regions[{{nt},{ne},{},{}}]
			mask = torch.squeeze(mask)

			mask_indxs = torch.nonzero(mask)
			r = mask_indxs[{{},{1}}]
			c = mask_indxs[{{},{2}}]

			centers[{{1},{nt},{ne}}] = torch.mean(r:double())
			centers[{{2},{nt},{ne}}] = torch.mean(c:double())

			areas[{{nt},{ne}}] = r:size(1)

		end
	end
end

-- Define a filters class here.

-- And then do something like:
-- filters.regions = regions
-- filters.centers = centers
-- filters.areas = areas

Filters = { regions = torch.Tensor(N_theta,N_e,visual_field_width,visual_field_width),
centers = torch.Tensor(2,N_theta,N_e),
areas = torch.Tensor(N_theta,N_e),
offsets_r = {},
offsets_c = {},
weights = {},
uniq_pix = {},
}
-- Don't allocate anything here


filters = Filters
nt=1
ne=1
for nt=1,N_theta do
	filters.offsets_r[nt] = {}
	filters.offsets_c[nt] = {}
	filters.weights[nt] = {}
	filters.uniq_pix[nt] = {}
	for ne=1,N_e do
		filters.offsets_r[nt][ne] = torch.Tensor()
		filters.offsets_c[nt][ne] = torch.Tensor()
		filters.weights[nt][ne] = torch.Tensor()
		filters.uniq_pix[nt][ne] = torch.Tensor()
		--print(nt,ne)
	end
end

filters.regions = regions
filters.center = centers
filters.areas = areas

-- We want to preserve smooth curves
blindspot_threshold = 0.0;

-- Initialize filters_offsets, filter weights, and 
-- Double check to see if this is working
center_r = math.floor(0.5+filters.regions:size(3)/2)
center_c = math.floor(0.5+filters.regions:size(4)/2)

for nt=1,N_theta do
	for ne=1,N_e do
		rc_spot = torch.nonzero(torch.squeeze(filters.regions[{{nt},{ne},{},{}}]))
		if rc_spot:nDimension()~=0 then
			--print(rc_spot:size())
			-- Redefine the Tensor dimensions:
			filters.offsets_r[nt][ne] = torch.Tensor(rc_spot:size(1))
			filters.offsets_c[nt][ne] = torch.Tensor(rc_spot:size(1))
			filters.weights[nt][ne] = torch.Tensor(rc_spot:size(1))
			filters.uniq_pix[nt][ne] = torch.Tensor(rc_spot:size())
			-- Assign new Tensor values:
			filters.offsets_r[nt][ne] = rc_spot[{{},{1}}] - center_r
			filters.offsets_c[nt][ne] = rc_spot[{{},{2}}] - center_c
			--idxs = filters.regions[]
			--filters.weights[nt][ne] = filters.regions[{{nt},{ne},{rc_spot[{{},{1}}]},{rc_spot[{{},{2}}]}}]
			-- Get Uniq pixels for the pooling region:
			--filters.uniq_pix[nt][ne][1] = rc_spot[{{},{1}}]
			--filters.uniq_pix[nt][ne][2] = rc_spot[{{},{2}}]
			filters.uniq_pix[nt][ne] = rc_spot
		end
	end
end

-- Optional Create Color palette: 
-- Updating this later:


-- Foveal Trimming:
-- Here we discard cells withing the fovea

fovea_radius = fovea
peripheral_filters = filters
n_e_limit = 1

for n_e=1,N_e do
	offsets_r = torch.abs(peripheral_filters.offsets_r[1][n_e])
	offsets_c = torch.abs(peripheral_filters.offsets_c[1][n_e])

	-- See how much of this cell is within the fovea:
	within = torch.le(offsets_r,fovea_radius/deg_per_pixel) and torch.le(offsets_c,fovea_radius/deg_per_pixel)

	if torch.sum(within)/within:size(1)<.5 then
		n_e_limit = n_e
		break
	end
end

-- 

-- Assign the new values:
peripheral_filters.regions = peripheral_filters.regions[{{},{n_e_limit,N_e},{},{}}]
peripheral_filters.centers = peripheral_filters.centers[{{},{},{n_e_limit,N_e}}]
peripheral_filters.areas = peripheral_filters.areas[{{},{n_e_limit,N_e}}]

-- Work on this

for n_e2=1,n_e_limit-1 do
	for n_t=1,N_theta do
		table.remove(peripheral_filters.offsets_r[n_t],1)
		table.remove(peripheral_filters.offsets_c[n_t],1)
		table.remove(peripheral_filters.weights[n_t],1)
		table.remove(peripheral_filters.uniq_pix[n_t],1)
	end
end
--peripheral_filters.offsets_r = peripheral_filters.offsets_r[{{},{n_e_limit,N_e},{},{}}]
--peripheral_filters.offsets_c = peripheral_filters.offsets_c[{{},{n_e_limit,N_e},{},{}}]
--peripheral_filters.weights = peripheral_filters.weights[{{},{n_e_limit,N_e},{},{}}]
--peripheral_filters.uniq_pix = peripheral_filters.uniq_pix[{{},{n_e_limit,N_e},{},{}}]
