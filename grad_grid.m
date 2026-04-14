function [ Gn_h , Ge_h , Gn_w , Ge_w , GRAD_grid_file ] = grad_grid ( indir_GRAD_grid, GRAD_grid_file , mjd , lat , lon , grid_res )
%
% grad_grid.m
%
% This routine bilinearly interpolates the GRAD gradients from the
% surrounding grid points to the desired site.
% http://vmf.geo.tuwien.ac.at/trop_products/GRID/
%
% On the temporal scale, the values from the two surrounding NWM epochs are
% linearly interpolated to the respective mjd. In the horizontal, a bilinear 
% interpolation is done.
%
% Reference for GRAD:
% Landskron, D. & Bohem, J. J Geod (2018). https://doi.org/10.1007/s00190-018-1127-1
%
%
% INPUT:
%         o indir_GRAD_grid ... input directory where the grid-wise GRAD files are stored
%         o GRAD_grid_file .... cell containing filenames and GRAD data, which is always passed with the function; must be set to '[]' by the user in the initial run
%         o mjd ............... modified Julian date
%         o lat ............... ellipsoidal latitude in radians
%         o lon ............... ellipsoidal longitude in radians
%         o grid_res........... grid resolution in degrees (possible: 1 or 5)
%
% OUTPUT:
%         o Gn_h .............. hydrostatic north gradient
%         o Ge_h .............. hydrostatic east gradient
%         o Gn_w .............. wet north gradient
%         o Ge_w .............. wet east gradient
%         o GRAD_grid_file .... cell containing filenames and GRAD data, which is always passed with the function; must be set to '[]' by the user in the initial run
%
% -------------------------------------------------------------------------
%
% written by Daniel Landskron (2018/02/21)
%
% =========================================================================



% save lat and lon also in degrees
lat_deg = lat*180/pi;
lon_deg = lon*180/pi;

% due to numerical issues, it might happen that the above conversion does not give exact results, e.g. in case of rad2deg(deg2rad(60)); in order to prevent this, lat_deg and lon_deg are rounded to the 10th decimal place
lat_deg = round(lat_deg,10);
lon_deg = round(lon_deg,10);


%% (1) convert the mjd to year, month, day in order to find the correct files


% find the two surrounding epochs
if mod(mjd,0.25)==0
    mjd_all = mjd;
else
    mjd_int = floor(mjd*4)/4 : 0.25 : ceil(mjd*4)/4;
    mjd_all = [mjd mjd_int];
end


hour = floor((mjd_all-floor(mjd_all))*24);   % get hours
minu = floor((((mjd_all-floor(mjd_all))*24)-hour)*60);   % get minutes
sec = (((((mjd_all-floor(mjd_all))*24)-hour)*60)-minu)*60;   % get seconds

% change secs, min hour whose sec==60
minu(sec==60) = minu(sec==60)+1;
hour(minu==60) = hour(minu==60)+1;
mjd_all(hour==24)=mjd_all(hour==24)+1;

% calc jd (yet wrong for hour==24)
jd_all = mjd_all+2400000.5;

% integer Julian date
jd_all_int = floor(jd_all+0.5);

aa = jd_all_int+32044;
bb = floor((4*aa+3)/146097);
cc = aa-floor((bb*146097)/4);
dd = floor((4*cc+3)/1461);
ee = cc-floor((1461*dd)/4);
mm = floor((5*ee+2)/153);

day = ee-floor((153*mm+2)/5)+1;
month = mm+3-12*floor(mm/10);
year = bb*100+dd-4800+floor(mm/10);

epoch = (mjd_all-floor(mjd_all))*24;

% derive related GRAD filename(s)
if length(mjd_all)==1   % if the observation epoch coincides with an NWM epoch
    filename = ['GRAD_' num2str(year(1)) sprintf('%02s',num2str(month(1))) sprintf('%02s',num2str(day(1))) '.H' sprintf('%02s',num2str(epoch(1)))];
else
    for i_mjd = 2:length(mjd_all)
        filename(i_mjd-1,:) = ['GRAD_' num2str(year(i_mjd)) sprintf('%02s',num2str(month(i_mjd))) sprintf('%02s',num2str(day(i_mjd))) '.H' sprintf('%02s',num2str(epoch(i_mjd)))];
    end
end

% only positive longitude in degrees
if lon_deg < 0
    lon = lon + 2*pi;
    lon_deg = (lon_deg + 360);
end



%% (2) check if new files have to be loaded or if the overtaken ones are sufficient


if isempty(GRAD_grid_file)   % in the first run, 'GRAD_file' is always empty and the orography_ell file has to loaded
    load_new = 1;
    GRAD_grid_file{1} = filename;   % replace the empty cell by the current filenames
    GRAD_grid_file{3} = indir_GRAD_grid;   % replace the empty cell by the current indir_GRAD_grid
    GRAD_grid_file{4} = [lat lon];
elseif strcmpi(GRAD_grid_file{1},filename)   &&   (lat > GRAD_grid_file{4}(1)   ||   (lat == GRAD_grid_file{4}(1) && lon <= GRAD_grid_file{4}(2) && lon >= grid_res/2))   &&   strcmpi(indir_GRAD_grid,GRAD_grid_file{3})   % if the current filenames are the same as in the forwarded files, and the coordinates are the same as well   
    load_new = 0;
    GRAD_data_all = GRAD_grid_file{2};
else   % if new files are required, then everything must be loaded anew
    load_new = 1;
    GRAD_grid_file{1} = filename;
    GRAD_grid_file{3} = indir_GRAD_grid;
    GRAD_grid_file{4} = [lat lon];
end



%% (3) find the indices of the 4 surrounding grid points


% find the coordinates (lat,lon) of the surrounding grid points
lat_all = 90-grid_res/2 : -grid_res : -90;
lon_all = 0+grid_res/2 : grid_res : 360;

% find the 2 closest latitudes
lat_temp = lat_deg-lat_all;
[~,ind_lat_int(1)] = min(abs(lat_temp));
ind_lat_int(2) = ind_lat_int(1)-sign(lat_temp(ind_lat_int(1)));

% find the two closest longitudes
lon_temp = lon_deg-lon_all;
[~,ind_lon_int(1)] = min(abs(lon_temp));
ind_lon_int(2) = ind_lon_int(1)+sign(lon_temp(ind_lon_int(1)));

% correct indices out of range
for i_ind = 1:2
    if ind_lat_int(i_ind)>length(lat_all); ind_lat_int(i_ind) = length(lat_all);                    end
    if ind_lat_int(i_ind)<1;               ind_lat_int(i_ind) = 1;                                  end
    if ind_lon_int(i_ind)>length(lon_all); ind_lon_int(i_ind) = ind_lon_int(i_ind)-length(lon_all); end
    if ind_lon_int(i_ind)<1;               ind_lon_int(i_ind) = ind_lon_int(i_ind)+length(lon_all); end
end

% define the indices
index(1) = (ind_lat_int(1)-1)*length(lon_all)+ind_lon_int(1);
index(2) = (ind_lat_int(1)-1)*length(lon_all)+ind_lon_int(2);
index(3) = (ind_lat_int(2)-1)*length(lon_all)+ind_lon_int(1);
index(4) = (ind_lat_int(2)-1)*length(lon_all)+ind_lon_int(2);



%% (4) read the correct data and perform a linear time interpolation from the surrounding two epochs
% read in with textscan, but only up to maximum index, everything before will be treated as headerlines


if load_new == 1
    
    for i_file = 1:size(filename,1)
        
        % read the files and collect the data
        if length(mjd_all)==1   % if the observation epoch coincides with an NWM epoch
            dat = fopen([indir_GRAD_grid '/' num2str(year(1)) '/' filename(1,:)]);
        else
            dat = fopen([indir_GRAD_grid '/' num2str(year(i_file+1)) '/' filename(i_file,:)]);
        end
        
        GRAD_data_all(i_file) = textscan(dat,'%f%f%f%f%f%f',max(index),'CommentStyle','!','CollectOutput',1);   % only read data up to the maximum index in order to save time
        fclose(dat);
        GRAD_grid_file{2} = GRAD_data_all;   % write the GRAD data to the forwarded variable
        GRAD_data{i_file} = GRAD_data_all{i_file}(index,:);   % reduce to the indices of the surrounding grid points
        
    end
else
    
    GRAD_data = cellfun(@(c) c(index,:),GRAD_data_all,'UniformOutput',false);   % reduce to the indices of the surrounding grid points
    
end


% initialize
GRAD_data_int = zeros(4,6);

% do the linear time interpolation for each argument; the results are the GRAD values for the surrounding grid points at the time of the measurement
iv_ind = 1:4;
if length(mjd_all)==1   % if the observation epoch coincides with an NWM epoch
    GRAD_data_int(iv_ind,1:6) = GRAD_data{1}(iv_ind,1:6);
else   % else perform the linear interpolation
    iv_line = 1:6;
    GRAD_data_int(iv_ind,iv_line) = GRAD_data{1}(iv_ind,iv_line) + (GRAD_data{2}(iv_ind,iv_line)-GRAD_data{1}(iv_ind,iv_line))*(mjd-mjd_int(1))/(mjd_int(2)-mjd_int(1));   % the appendix 'h0' means that the values are valid at zero height
end




%% (5) perform the bilinear interpolation


if length(unique(index)) == 1   % if the point is directly on a grid point
    
    Gn_h = GRAD_data_int(1,3);
    Ge_h = GRAD_data_int(1,4);
    Gn_w = GRAD_data_int(1,5);
    Ge_w = GRAD_data_int(1,6);
    
else
    
    % bilinear interpolation (interpreted as two 1D linear interpolations for lat and lon, but programmed without subfunctions)

    % (a) linear interpolation for longitude
    if ~isequal(GRAD_data_int(1,2), GRAD_data_int(2,2))   % if longitude must be interpolated
        Gn_h_lon1 = GRAD_data_int(1,3) + (GRAD_data_int(2,3)-GRAD_data_int(1,3))*(lon_deg-GRAD_data_int(1,2))/(GRAD_data_int(2,2)-GRAD_data_int(1,2));
        Gn_h_lon2 = GRAD_data_int(3,3) + (GRAD_data_int(4,3)-GRAD_data_int(3,3))*(lon_deg-GRAD_data_int(3,2))/(GRAD_data_int(4,2)-GRAD_data_int(3,2));
        Ge_h_lon1 = GRAD_data_int(1,4) + (GRAD_data_int(2,4)-GRAD_data_int(1,4))*(lon_deg-GRAD_data_int(1,2))/(GRAD_data_int(2,2)-GRAD_data_int(1,2));
        Ge_h_lon2 = GRAD_data_int(3,4) + (GRAD_data_int(4,4)-GRAD_data_int(3,4))*(lon_deg-GRAD_data_int(3,2))/(GRAD_data_int(4,2)-GRAD_data_int(3,2));
        Gn_w_lon1 = GRAD_data_int(1,5) + (GRAD_data_int(2,5)-GRAD_data_int(1,5))*(lon_deg-GRAD_data_int(1,2))/(GRAD_data_int(2,2)-GRAD_data_int(1,2));
        Gn_w_lon2 = GRAD_data_int(3,5) + (GRAD_data_int(4,5)-GRAD_data_int(3,5))*(lon_deg-GRAD_data_int(3,2))/(GRAD_data_int(4,2)-GRAD_data_int(3,2));
        Ge_w_lon1 = GRAD_data_int(1,6) + (GRAD_data_int(2,6)-GRAD_data_int(1,6))*(lon_deg-GRAD_data_int(1,2))/(GRAD_data_int(2,2)-GRAD_data_int(1,2));
        Ge_w_lon2 = GRAD_data_int(3,6) + (GRAD_data_int(4,6)-GRAD_data_int(3,6))*(lon_deg-GRAD_data_int(3,2))/(GRAD_data_int(4,2)-GRAD_data_int(3,2));
    else   % if the station coincides with the longitude of the grid
        Gn_h_lon1 = GRAD_data_int(1,3);
        Gn_h_lon2 = GRAD_data_int(3,3);
        Ge_h_lon1 = GRAD_data_int(1,4);
        Ge_h_lon2 = GRAD_data_int(3,4);
        Gn_w_lon1 = GRAD_data_int(1,5);
        Gn_w_lon2 = GRAD_data_int(3,5);
        Ge_w_lon1 = GRAD_data_int(1,6);
        Ge_w_lon2 = GRAD_data_int(3,6);
    end
    
    % linear interpolation for latitude
    if ~isequal(GRAD_data_int(1,1), GRAD_data_int(3,1))   % if latitude must be interpolated
        Gn_h = Gn_h_lon1 + (Gn_h_lon2-Gn_h_lon1)*(lat_deg-GRAD_data_int(1,1))/(GRAD_data_int(3,1)-GRAD_data_int(1,1));
        Ge_h = Ge_h_lon1 + (Ge_h_lon2-Ge_h_lon1)*(lat_deg-GRAD_data_int(1,1))/(GRAD_data_int(3,1)-GRAD_data_int(1,1));
        Gn_w = Gn_w_lon1 + (Gn_w_lon2-Gn_w_lon1)*(lat_deg-GRAD_data_int(1,1))/(GRAD_data_int(3,1)-GRAD_data_int(1,1));
        Ge_w = Ge_w_lon1 + (Ge_w_lon2-Ge_w_lon1)*(lat_deg-GRAD_data_int(1,1))/(GRAD_data_int(3,1)-GRAD_data_int(1,1));
    else   % if the station coincides with the latitude of the grid
        Gn_h = Gn_h_lon1;
        Ge_h = Ge_h_lon1;
        Gn_w = Gn_w_lon1;
        Ge_w = Ge_w_lon1;
    end
      
end