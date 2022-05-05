clear;clc;

% anomoly
size = 10;
anomaly = gausswin(size) * gausswin(size)';
element = [anomaly, -anomaly; -anomaly, anomaly];
vpmodel = repmat(element, 5, 20);
vpmodel = vpmodel/ max(vpmodel(:));

% velocity model
vp_background = ones(101, 401) * 3500;
vp_checkboard = vp_background;
vp_checkboard(1:end-1, 1:end-1) = vp_checkboard(1:end-1, 1:end-1) + vpmodel * 3500 * 0.10;


vp_checkboard = vp_checkboard';
vp_background = vp_background';

save background_401_201_25m.dat vp_background -ascii
save checkboard_801_201_25m.dat vp_checkboard -ascii