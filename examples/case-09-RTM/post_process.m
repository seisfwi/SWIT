clear;clc;close all
addpath(genpath('/home/haipeng/Nutstore Files/Nutstore/Packages/SWIT-process/'), '-begin');

nx = 1201;
nz = 515;
dx = 10;
z= [0:nz-1]*dx;

% mask
vp = load('./model/SEAM_I_1201_515_10m.dat')';
mask = vp==1490;

% image
fp = fopen('./outputs/gradient/RTM-image.bin', 'rb'); RTM = fread(fp, 'float32'); fclose(fp);
RTM = reshape(RTM, [nz, nx]);
RTM = RTM/max(abs(RTM(:)));

% forward illum
fp = fopen('./data/syn/src0_illum_forw.bin', 'rb'); illum = fread(fp, 'float32'); fclose(fp);
illum = reshape(illum, [nz, nx]);

% illum
RTM = RTM./ power(illum, 2.);
RTM(isnan(RTM)) = 0.0;

% bandpass
% RTM = butterband(RTM, z, 0.002, 0, 4, 0);

% diff
order = 2;
if order > 0
    RTM = diff(RTM,order);
    RTM = [zeros(order, nx); RTM];
end 

% mask 
RTM(mask) = 0.0;

% normlize
RTM = RTM / max(abs(RTM(:)));

figure(1)
imagesc(RTM); colormap(gray)
caxis([-1 1] * 0.1)


save RTM_515_1201.dat RTM -ascii 
% print(gcf, 'RTM-image-no-diff.png', '-dpng', '-r300');

