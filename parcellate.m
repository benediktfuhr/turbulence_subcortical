function [ts] = parcellate(imgfile, parcfile)

%% Function that takes as an input the rest-State file and the resampled Schaefer
%  File and computes the Schaefer parc

%Read in data
rest_data = read_avw(imgfile);
atlas = read_avw(parcfile);

%% Time-Series Extraction
region_labels = unique(atlas(:));
region_labels(region_labels == 0) = [];
n_regions = length(region_labels);
ntp = size(rest_data,4);

% tic; 

% pre-allocate
ts = nan(n_regions, ntp);

% loop over regions
for iParcel = 1:n_regions

    mask = atlas == iParcel; % Create the parcel mask

    % mike's method -----------------------------------------------------------
    % loop over timepoints
    for itp = 1:ntp
        tmp_data = rest_data(:,:,:,itp);
        parcel_voxels = tmp_data(mask);
        ts(iParcel,itp) = nanmean(parcel_voxels);
    end % for itp
    % -------------------------------------------------------------------------

    % % benedikt's method -----------------------------------------------------
    % parcel_voxels = rest_data(repmat(mask, [1 1 1 ntp])); % Extract all parcel voxels for all time points
    % ts(iParcel, :) = nanmean(reshape(parcel_voxels, [], ntp), 1); % Average over spatial dimensions
    % -------------------------------------------------------------------------
    
    % Normalize time series
    ts(iParcel, :) = (ts(iParcel, :) - nanmean(ts(iParcel, :))) / std(ts(iParcel, :));
    
end % for iParcel
% toc;