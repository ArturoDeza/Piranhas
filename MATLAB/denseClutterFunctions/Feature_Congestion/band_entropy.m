function en_band = band_entropy(map, wlevels, wor)

% en_band = band_entropy(map, wlevels, wor)
% Inputs: 
%     "map": a monochromatic image
%     "wlevels": the number of spatial scales for the subband decomposition
%     "wor": the number of orientations for the subband decomposition
% Outputs:
%     "en_band": a vector containing Shannon entropies of all the subbands

% Reference: Ruth Rosenholtz, Yuanzhen Li, and Lisa Nakano. "Measuring
% Visual Clutter". To appear in Journal of Vision, 7(2).

% Ruth Rosenholtz, Yuanzhen Li, and Lisa Nakano, March 2007.

% decompose the image into subbands:
[C,S] = buildSFpyr(map,wlevels,wor-1);
ind = prod(S(1,:));
band_ind = 2;
for ii = 1:wlevels
    for jj = 1:wor
        band = C(ind+1:ind+prod(S(band_ind,:)));
        ind = ind+prod(S(band_ind,:));
        % Shannon entropy:
        en_band(1,band_ind-1) = entropy(band(:)'); 
        band_ind = band_ind + 1;
    end    
end
band = C(ind+1:ind+prod(S(band_ind,:)));
ind = ind+prod(S(band_ind,:));
en_band(1,band_ind-1) = entropy(band(:)'); 
