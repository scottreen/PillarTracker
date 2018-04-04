% This script deomo to how to create an interpolated force map.
%% load data and metadata from binary file. change the file path if necessary
[meta, data] = read_pillar_bin('C:\Users\xiaochun\Desktop\MDA_pillar_10sec_15min__R3D.bin');
pixel_size = meta.pixel_size;

%% interpolation on specified frame
frame = 1;
dx = data.DX(frame, :);
dy = data.DY(frame, :);
d = sqrt(dx.^2 + dy.^2);
% grid x and y
x = data.tracksX(frame, :)-dx;
y = data.tracksY(frame, :)-dy;

% remove NaN
flags_nan = (isnan(x) | isnan(y)); 
x = x(~flags_nan);
y = y(~flags_nan);
d = d(~flags_nan);

% range of x and y
xmin = floor(min(x));
xmax = ceil(max(x));
ymin = floor(min(y));
ymax = ceil(max(y));

% interpolation
[xq, yq] = meshgrid(xmin:xmax, ymin:ymax);
dq = griddata(x, y, d, xq, yq, 'cubic');
dq(dq<0) = 0;

% compute the force map, for example, 500nm diameter, and 800nm high,
% 100kPa Young's Modulus. the force is calculated in pN.
D = 500;  %pillar diameter, in nm
L = 800;  %pillar tall or length, in nm
E = 0.1;  %Young's modulus, in pN/nm^2, 100kPa
F = 3*pi/64*D^4/L^3*E*dq*pixel_size; % in pN
% plot force map
figure, imagesc(F), colorbar
title('Force Map (pN)');

