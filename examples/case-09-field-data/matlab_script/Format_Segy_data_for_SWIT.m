%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Helpful tools for SWIT-1.0
% By Haipeng Li @ USTC
% Aug. 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc; close all;
addpath(genpath('./SegyMAT/'), '-begin');

% read your data
segyfile = './shot020.segy';
[data, TH, H] = ReadSegy(segyfile);


% set the desired receiver range (same with the FWI code)
rec_x_beg = 0.0;
rec_x_end = 50e3;
rec_dx    = 25;
recX_new  = [rec_x_beg : rec_dx : rec_x_end];

% get basic parameters from data
nt = size(data, 1);
dt = TH(1).dt * 1e-6;
srcX = cell2mat({TH.SourceX}); srcX = srcX(1);
recX = cell2mat({TH.GroupX});

% write a empty Segy file and then convert to SU file
WriteSegy('empty.segy', zeros(nt, length(recX_new)), 'dt', dt);
Segy2Su('empty.segy');

% load the empty SU file
[data_su, TH_su, H_su] = ReadSu('empty.segy.su');

% put the original data into this empty data
for itrace = 1 : length(recX)
    % find nearest receiver position
    [~, index] = min(power(recX_new - recX(itrace), 2));
    data_su(:, index) = data(:, itrace);
end
% add necessary headers
for itrace = 1 : length(recX_new)
    TH_su(itrace).SourceX = srcX;
    TH_su(itrace).GroupX = recX_new(itrace);
    % add more header here as you like 
end

% remove the dummy files
system('rm empty.segy empty.segy.su');

% write the final SU file
WriteSuStructure([segyfile(1:end-4), 'su'], H_su, TH_su, data_su);



