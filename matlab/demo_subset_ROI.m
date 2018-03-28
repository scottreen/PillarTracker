%% load csv file with specified ROI. change the file path if necessary
subset_table = readtable('C:\Users\xiaochun\Desktop\Data for Frame=1.csv');
PillarIndex=subset_table.PillarIndex;

%% load data and metadata from binary file. change the file path if necessary
[meta, data] = read_pillar_bin('I:\mingxi\casminus_1.3um_2dhex_cell03_R3D\casminus_1.3um_2dhex_cell03_R3D.bin');
pixel_size = meta.pixel_size;

%% get subset data
%remove the NaN index
subset = PillarIndex(~isnan(PillarIndex));
%tracks subset
sub_X = data.tracksX(:, subset);
sub_Y = data.tracksY(:, subset);
%deflections subset
sub_DX = data.DX(:, subset);
sub_DY = data.DY(:, subset);

%% compute the deflections, use the Euclidean distance
sub_DIS = sqrt(sub_DX.^2 + sub_DY.^2);

%% convert the pixel unit to nm by multiplying the pixel_size
sub_deflections = sub_DIS*pixel_size;

%% plot the first frame
scatter(data.tracksX(1,:), data.tracksY(1,:),'x');
hold on, scatter(sub_X(1,:), sub_Y(1,:),'x'); hold off

%% compute the force
D = 500;  %pillar diameter, in nm
L = 800;  %pillar tall or length, in nm
E = 0.1;  %Young's modulus, in pN/nm^2, 100kPa
% d = 35;    %deflection, in nm
% force calculation in pN
% F = 3*pi/64*D^4/L^3*E*d; % ~62.9 pN
F = 3*pi/64*D^4/L^3*E*sub_deflections;
