%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate horizontal gradient from grid-wise GRAD product with height correction
%Creat: Yaozong Zhou, GNSS Research Center, Wuhan University
%Date: 2026-03-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%----------------------------------USER INPUT------------------------------------
%
%start and end time, format: YYYYMMDDHH
symdh='2019080911';
eymdh='2019080911';
%
%time resolution, unit: hour
res_hour=1;
%
%GRAD product dir, automatic recognition year
grad_dir='./';
%
%coordinate file, format: station name, B (latitude), L (longitue), H (geodetic height)
blhfile='Example.BLH';
%
%result file, format: cmjd, gnh, gnh_hc, geh, geh_hc, geh, geh_hc, gew, gew_hc, unit: mm
resfile='Example_Result.HG';
%
%---------------------------------CONVERT TIME----------------------------------
%
sarr=sscanf(symdh,'%04d%02d%02d%02d');
earr=sscanf(eymdh,'%04d%02d%02d%02d');
sy=sarr(1); sm=sarr(2); sd=sarr(3);shr=sarr(4);
ey=earr(1); em=earr(2); ed=earr(3);eh=earr(4);
smjd=modified_julday(sd,sm,sy);
emjd=modified_julday(ed,em,ey);
%
%-------------------------------READ COORDINATE---------------------------------
%
[sta,sb,sl,sh] = read_blhfile(blhfile);
nsta = length(sta);
consrad=180/pi;
rsb = sb/consrad;
rsl = sl/consrad;
%
%-------------------------------MJD AND HOUR LOOP-------------------------------
%
for imjd=smjd:emjd
   %
   %Decide start and end hours
   ish=0;
   ieh=23;
   if (imjd==smjd)
   if (shr>0)
     ish=shr; 
   end
   end
   if (imjd==emjd)
   if (eh<23)
     ieh=eh; 
   end
   end
   %  
   %For every hour
   for ihour=ish:ieh
     if (rem(ihour,res_hour) == 0)
       % 
       %Get current mjd
       cmjd=imjd+ihour/24;
       %
       %--------------------------STATION LOOP----------------------------------
       %
       for ista=1:nsta
           %
           %Set initial value
           GRAD_grid_file = [];
           %
           %---------------EXTRACT HORIZONTAL GRADIENT-------------------------
           %
           %Product downloaded from https://vmf.geo.tuwien.ac.at/trop_products/GRID/1x1/GRAD/GRAD_OP/
           %Function grad_grid provided by TU WIEN in https://vmf.geo.tuwien.ac.at/codes/
           [ gnh , geh , gnw , gew , GRAD_grid_file ] = grad_grid ( grad_dir, GRAD_grid_file , cmjd , rsb(ista), rsl(ista) , 1 );
           %
           %-------------------HEIGHT CORRECTION-------------------------------
           %
           %FIRST STEP: calculate model parameters
           %
           %grid height: grid_hgh
           %first-order exponential coefficients: gnh_C1,geh_C1,gnw_C1,gew_C1
           [grid_hgh,gnh_C1,geh_C1,gnw_C1,gew_C1] = wtm_hc_hg (cmjd, rsb(ista), rsl(ista), 0);
           %
           %SECOND STEP: height correction application
           %
           %Calculate height difference between station and grid
           hk=(sh(ista)-grid_hgh)/1000;
           %Correction, when height difference larger than 1 km
           if (hk>=1)
             gnh_hc=gnh*exp(hk*gnh_C1);
             geh_hc=geh*exp(hk*geh_C1);
             gnw_hc=gnw*exp(hk*gnw_C1);
             gew_hc=gew*exp(hk*gew_C1);
           else
             gnh_hc=gnh;
             geh_hc=geh;
             gnw_hc=gnw;
             gew_hc=gew;            
           end
           %
           %-------------------------OUTPUT-------------------------------------
           %
           %Define the result files
           resf=fopen(resfile,'a');
	   	   %
	       %Write time and SPDs
           fprintf(resf,'%10.2f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n',cmjd,gnh,gnh_hc,geh,geh_hc,geh,geh_hc,gew,gew_hc);
           %
           %Close file
           fclose(resf);
       end
     end
   end
end