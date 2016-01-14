%% Load HTP Dataset and perform tracking 

% unsing the Crocker, Weeks, and Grier Algorithm (http://www.physics.emory.edu/rweeks/idl/index.html)

% clear, close all, clc, clear
%%  %%%%%%%%%% INPUT Parameters %%%%%%%%%%
function [res]=CG_tracking(filename_peaks);

% filename_peaks='A647_EGF_10ms_1500mW_COT_Au__2_MMStack_locResults_DC';     % filename of TS output file

max_disp = 14;  % in pxl
min_pos = 1;    % good - eliminate if fewer than good valid positions
gap = 2280;     % mem - number of time steps that a particle can be 'lost' and then recovered again
quiet = 1;      % quiet - 1 = no text

%% Load Data

filename_peaks2=[filename_peaks '.dat'];
peaks=dlmread(filename_peaks2,',',1,0);

% 
% file = fopen(filename_peaks2);
% line = fgetl(file);
% h = regexp( line, ',', 'split' );
% 
% x = strmatch('x [nm]',h);
% y = strmatch('y [nm]',h);
% frame = strmatch('frame',h);

% Create pos_list for track.m

pos_list(:,1)=peaks(:,1);                   % in pxl
pos_list(:,2)=peaks(:,2);                   % in pxl
pos_list(:,3)=peaks(:,4);               % dt in seconds

fprintf('\n -- Data Loaded --\n')

%% Track unsing the Crocker, Weeks, and Grier Algorithm (http://www.physics.emory.edu/rweeks/idl/index.html)

param=struct('mem',gap,'dim',2,'good',min_pos,'quiet',quiet);
res=trackGT(pos_list,max_disp,param); % variable XYT, maximum displacement in pxl

fprintf('\n -- Tracking Done --\n')

%% Plot All Tracks

% scatter(res(:,1),res(:,2),5,res(:,4));

% for i=1:max(res(:,4))
%     
%     vx=find(res(:,4)==i);
%     
%     plot(res(vx,1),res(vx,2));hold on;
% 
% end

%% Save tracks

filenamec1=['Tracks_' filename_peaks];

save(filenamec1,'res');

fprintf('\n -- Tracks Saved --\n')

% 1 - x
% 2 - y
% 3 - time in seconds
% 4 - track ID
end