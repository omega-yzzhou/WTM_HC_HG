function [mjd]=modified_julday(iday,imonth,iyear)
%
% modified_julday.m
%
% GNSS Research Center, Wuhan University, 2026 
%
% This subroutine calculates modified julian day according to year, month and day
%
% INPUT: 
%        o iyear: year 4-digits
%        o imonth: 1 to 12 
%        o iday: 1 to 31
% OUTPUT:
%        o mjd: modified julian day
%----------------------------------------------------------------------------
% 
% m-file created by Yaozong Zhou
%
%----------------------------------------------------------------------------
% History:
%
% 2026-03-14: m-file created
%
%----------------------------------------------------------------------------
%
%day of year for month
doy_of_month=[0,31,59,90,120,151,181,212,243,273,304,334];
%
%check date
if(iyear<0 || imonth<0 || iday<0 || imonth>12 || iday>366 || (imonth ~=0 && iday>31))
    str=['ERROR(modified_julday): incorrect date (year,month,day):',...
        num2str(iyear),' ',num2str(imonth),' ',num2str(iday)];
    disp(str);
    exit(1);
end
%
%suit for month smaller than 2
iyr=iyear;
if(imonth<=2) 
    iyr=iyr-1; 
end
%
%convert year and day to day
mjd=365*iyear-678941+fix(iyr/4)-fix(iyr/100)+fix(iyr/400)+iday;
%
%convert month to day
if(imonth~=0) 
    mjd=mjd+doy_of_month(imonth); 
end