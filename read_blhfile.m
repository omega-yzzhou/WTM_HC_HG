function [sta,sb,sl,sh] = read_blhfile(blhfile)
%
% read_blhfile.m
%
% GNSS Research Center of Wuhan University, 2020 
%
% This subroutine reads the .BLH file
%
% INPUT: 
%        o blhfile: .BLH file path and name
%
% OUTPUT:
%        o sta: station name
%        o sb: station latitude
%        o sl: station longitude
%        o sh: station height
%
%--------------------------------------------------------------------------
% 
% m-file created by Yaozong Zhou
%
%--------------------------------------------------------------------------
% History:
%
% 2026-03-14: m-file created
%
%--------------------------------------------------------------------------
%    
%Open .BLH file
inpf=fopen(blhfile,'r');
%
%Read .BLH file
inpc=textscan(inpf,'%s %f %f %f','headerlines',0);
sta=inpc{1};
sb=inpc{2};
sl=inpc{3};
sh=inpc{4};
%
%Close file
fclose(inpf);