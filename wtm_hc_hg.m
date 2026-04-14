function [grid_hgh,gnh_C1,geh_C1,gnw_C1,gew_C1] = wtm_hc_hg (mjd,lat,lon,it)
%
% wtm_hc_hg.m
%
% GNSS Research Center, Wuhan University, 2026 
%
% This subroutine calculates grid height as well as first-order exponential 
% coefficients from Wuhan University Tropospheric Model (WTM) for horizontal
% gradient height correction.
%
% INPUT: 
%        o mjd  - modified Julian Date (scalar)
%        o lat  - ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
%        o lon  - longitude in radians [-pi:pi] or [0:2pi] (vector)
%        o it   - case 0: with annual/semiannual variation; case 1: static quantities
%
% OUTPUT:
%        o grid_hgh - grid height at the station location (vector) (unit: m)
%        o gnh_C1   - hydrostatic north gradient coefficient (vector)
%        o geh_C1   - hydrostatic east gradient coefficient (vector)
%        o gnw_C1   - wet north gradient coefficient (vector)
%        o gew_C1   - wet east gradient coefficient (vector)
%----------------------------------------------------------------------------
% 
% m-file created by Yaozong Zhou
%
%----------------------------------------------------------------------------
% History:
%
% 2026-03-14: m-file created according to gpt3_1.m provided by TU WIEN in 
%             https://vmf.geo.tuwien.ac.at/codes/
%
%----------------------------------------------------------------------------
%
% read .MOD file
fid = fopen('WTM_HC_HG.MOD','r');
C = textscan( fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 1 , 'CollectOutput', true );
C = C{1};
fclose (fid);
%
% read grid height and coefficient
hgh_grid  = C(:,3);         % height
gnhc1_grid  = C(:,4:8);     % Gnh_C1
gehc1_grid  = C(:,9:13);    % Geh_C1
gnwc1_grid  = C(:,14:18);   % Gnw_C1
gewc1_grid  = C(:,19:23);   % Gew_C1
%
% convert mjd to doy
hour = floor((mjd-floor(mjd))*24);   % get hours
minu = floor((((mjd-floor(mjd))*24)-hour)*60);   % get minutes
sec = (((((mjd-floor(mjd))*24)-hour)*60)-minu)*60;   % get seconds
%
% change secs, min hour whose sec==60
minu(sec==60) = minu(sec==60)+1;
sec(sec==60) = 0;
hour(minu==60) = hour(minu==60)+1;
minu(minu==60)=0;
%
% calc jd (yet wrong for hour==24)
jd = mjd+2400000.5;
%
% if hr==24, correct jd and set hour==0
jd(hour==24)=jd(hour==24)+1;
hour(hour==24)=0;
%
% integer julian date
jd_int = floor(jd+0.5);
%
aa = jd_int+32044;
bb = floor((4*aa+3)/146097);
cc = aa-floor((bb*146097)/4);
dd = floor((4*cc+3)/1461);
ee = cc-floor((1461*dd)/4);
mm = floor((5*ee+2)/153);
%
day = ee-floor((153*mm+2)/5)+1;
month = mm+3-12*floor(mm/10);
year = bb*100+dd-4800+floor(mm/10);
%
% first check if the specified year is leap year or not (logical output)
leapYear = ((mod(year,4) == 0 & mod(year,100) ~= 0) | mod(year,400) == 0);
%
days = [31 28 31 30 31 30 31 31 30 31 30 31];
doy = sum(days(1:month-1)) + day;
if leapYear == 1 && month > 2
    doy = doy + 1;
end
doy = doy + mjd-floor(mjd);   % add decimal places
%
% factors for amplitudes
if (it==1) % then  constant parameters
    cosfy = 0;
    coshy = 0;
    sinfy = 0;
    sinhy = 0;
else 
    cosfy = cos(doy/365.25*2*pi);   % coefficient for A1
    coshy = cos(doy/365.25*4*pi);   % coefficient for B1
    sinfy = sin(doy/365.25*2*pi);   % coefficient for A2
    sinhy = sin(doy/365.25*4*pi);   % coefficient for B2
end
%
% determine the number of stations
nstat = length(lat);
%
% initialization
grid_hgh = zeros([nstat , 1]);
gnh_C1 = zeros([nstat , 1]);
geh_C1 = zeros([nstat , 1]);
gnw_C1 = zeros([nstat , 1]);
gew_C1 = zeros([nstat , 1]);
%
% loop over stations
for k = 1:nstat
    %
    % only positive longitude in degrees
    if lon(k) < 0
        plon = (lon(k) + 2*pi)*180/pi;
    else
        plon = lon(k)*180/pi;
    end
    % transform to polar distance in degrees
    ppod = (-lat(k) + pi/2)*180/pi; 
    %
    % find the index (line in the grid file) of the nearest point
    % changed for the 1 degree grid
    ipod = floor(ppod+1); 
    ilon = floor(plon+1);
    %
    % normalized (to one) differences, can be positive or negative
    % changed for the 1 degree grid
    diffpod = (ppod - (ipod - 0.5));
    difflon = (plon - (ilon - 0.5));
	% changed for the 1 degree grid
    if ipod == 181
        ipod = 180;
    end
    if ilon == 361
		ilon = 1;
    end
    if ilon == 0
        ilon = 360;
    end
    %
    % get the number of the corresponding line
    % changed for the 1 degree grid
    indx(1) = (ipod - 1)*360 + ilon;
    %
    % near the poles: nearest neighbour interpolation, otherwise: bilinear
    % with the 1 degree grid the limits are lower and upper
    bilinear = 0;
    if ppod > 0.5 && ppod < 179.5 
           bilinear = 1;          
    end           
    %
    % case of nearest neighbourhood
    if bilinear == 0
        %
        ix = indx(1);
        %    
        % calculate grid height and coefficients
        grid_hgh(k) = hgh_grid(ix,1);
        gnh_C1(k) = gnhc1_grid(ix,1) + gnhc1_grid(ix,2)*cosfy + gnhc1_grid(ix,3)*sinfy + gnhc1_grid(ix,4)*coshy + gnhc1_grid(ix,5)*sinhy;
        geh_C1(k) = gehc1_grid(ix,1) + gehc1_grid(ix,2)*cosfy + gehc1_grid(ix,3)*sinfy + gehc1_grid(ix,4)*coshy + gehc1_grid(ix,5)*sinhy;
        gnw_C1(k) = gnwc1_grid(ix,1) + gnwc1_grid(ix,2)*cosfy + gnwc1_grid(ix,3)*sinfy + gnwc1_grid(ix,4)*coshy + gnwc1_grid(ix,5)*sinhy;
        gew_C1(k) = gewc1_grid(ix,1) + gewc1_grid(ix,2)*cosfy + gewc1_grid(ix,3)*sinfy + gewc1_grid(ix,4)*coshy + gewc1_grid(ix,5)*sinhy;
		%
     else % bilinear interpolation
        %
        ipod1 = ipod + sign(diffpod);
        ilon1 = ilon + sign(difflon);
		% changed for the 1 degree grid
        if ilon1 == 361
            ilon1 = 1;
        end
        if ilon1 == 0
            ilon1 = 360;
        end
        %
        % get the number of the line
		% changed for the 1 degree grid
        indx(2) = (ipod1 - 1)*360 + ilon;  % along same longitude
        indx(3) = (ipod  - 1)*360 + ilon1; % along same polar distance
        indx(4) = (ipod1 - 1)*360 + ilon1; % diagonal
        %
        % calculate height and coefficients near four grid point (gp) 
        hgh_gp4 = hgh_grid(indx);
        gnhc1_gp4 = gnhc1_grid(indx,1) + gnhc1_grid(indx,2)*cosfy + gnhc1_grid(indx,3)*sinfy + gnhc1_grid(indx,4)*coshy + gnhc1_grid(indx,5)*sinhy;
        gehc1_gp4 = gehc1_grid(indx,1) + gehc1_grid(indx,2)*cosfy + gehc1_grid(indx,3)*sinfy + gehc1_grid(indx,4)*coshy + gehc1_grid(indx,5)*sinhy;
        gnwc1_gp4 = gnwc1_grid(indx,1) + gnwc1_grid(indx,2)*cosfy + gnwc1_grid(indx,3)*sinfy + gnwc1_grid(indx,4)*coshy + gnwc1_grid(indx,5)*sinhy;
        gewc1_gp4 = gewc1_grid(indx,1) + gewc1_grid(indx,2)*cosfy + gewc1_grid(indx,3)*sinfy + gewc1_grid(indx,4)*coshy + gewc1_grid(indx,5)*sinhy;	
        %    
        dnpod1 = abs(diffpod); % distance nearer point
        dnpod2 = 1 - dnpod1;   % distance to distant point
        dnlon1 = abs(difflon);
        dnlon2 = 1 - dnlon1;
		%
		% interpolation
        %grid_hgh
        R1 = dnpod2*hgh_gp4(1)+dnpod1*hgh_gp4(2);
        R2 = dnpod2*hgh_gp4(3)+dnpod1*hgh_gp4(4);
        grid_hgh(k) = dnlon2*R1+dnlon1*R2;
        %gnh_C1
        R1 = dnpod2*gnhc1_gp4(1)+dnpod1*gnhc1_gp4(2);
        R2 = dnpod2*gnhc1_gp4(3)+dnpod1*gnhc1_gp4(4);
        gnh_C1(k) = dnlon2*R1+dnlon1*R2;
        %geh_C1
        R1 = dnpod2*gehc1_gp4(1)+dnpod1*gehc1_gp4(2);
        R2 = dnpod2*gehc1_gp4(3)+dnpod1*gehc1_gp4(4);
        geh_C1(k) = dnlon2*R1+dnlon1*R2;
        %gnw_C1
        R1 = dnpod2*gnwc1_gp4(1)+dnpod1*gnwc1_gp4(2);
        R2 = dnpod2*gnwc1_gp4(3)+dnpod1*gnwc1_gp4(4);
        gnw_C1(k) = dnlon2*R1+dnlon1*R2;
        %gew_C1
        R1 = dnpod2*gewc1_gp4(1)+dnpod1*gewc1_gp4(2);
        R2 = dnpod2*gewc1_gp4(3)+dnpod1*gewc1_gp4(4);
        gew_C1(k) = dnlon2*R1+dnlon1*R2;           
    end 
end