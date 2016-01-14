clear, close all, clc, clear

%%  %%%%%%%%%%INPUT%%%%%%%%%%

filename_peaks='A647_EGF_10ms_1500mW_COT_Au__1_MMStack_locResults_DC';     % filename of TS output file
filename_peaks2=[filename_peaks '.dat'];
peaks=dlmread(filename_peaks2,',',1,0);


file = fopen(filename_peaks2);
line = fgetl(file);
h = regexp( line, ',', 'split' );

x = strmatch('x [nm]',h);
y = strmatch('y [nm]',h);
frame = strmatch('frame',h);

% Create pos_list for track.m

pos_list(:,1)=peaks(:,x)/100;                   % in pxl
pos_list(:,2)=peaks(:,y)/100;                   % in pxl
pos_list(:,3)=peaks(:,frame)*0.01;              % dt in seconds

fprintf('\n -- Data loaded --\n')

%% Track unsing the Crocker, Weeks, and Grier Algorithm (http://www.physics.emory.edu/rweeks/idl/index.html)

% mem - number of time steps that a particle can be 'lost' and then recovered again
% dim - default 2
% good - eliminate if fewer than good valid positions
% quiet - 1 = no text

param=struct('mem',10,'dim',2,'good',1,'quiet',1);
res=trackGT(pos_list,0.2,param); % variable XYT, maximum displacement in pxl

fprintf('\n -- Tracking done --\n')

%% Save tracks

filenamec1=['Tracks_' filename_peaks];

save(filenamec1,'res');

fprintf('\n -- Tracks Saved --\n')

% 1 - x
% 2 - y
% 3 - time in seconds
% 4 - track ID

%% Plot tracks longer than min_length
min_length=1;
max_length=10000;
image_name='STD_A647_EGF_10ms_1500mW_COT_Au__1_MMStack_Pos0.ome.tif';

plot_tracks(res,min_length,max_length,image_name,1); 

% tracks output, 
% minimum track length
% name of a file to overlay
% 1 = overlay; 0 = dont overlay

%% Track length histogram

tracklength=[];

for index=1:max(res(:,4))
    
    track=find(res(:,4)==index);
    tracklength=cat(1,tracklength,length(track));
    
    clear track
end

figure
hist(tracklength(tracklength>14),30)
% MeanTrackLength = mean(tracklength(tracklength>14))
mean(tracklength)
%% Generate variable tracks to interact with @msdanalyzer

min_length=5;       % min trajectory length
count=1;            % counter to select field
dx=0.107;           % conversion from pxl to mum

%%%%%%%%%%%%%%%%%%%%%% ROI selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%

upperx= max(res(:,1));         
lowerx=0;

uppery= max(res(:,2));
lowery= 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vx=find(res(:,1) < upperx & res(:,1) > lowerx);
subset=res(vx);
subset(:,2)=res(vx,2);
subset(:,3)=res(vx,3);
subset(:,4)=res(vx,4);

vy=find(subset(:,2) < uppery & subset(:,2) > lowery);
res_wROI=subset(vy);
res_wROI(:,2)=subset(vy,2);
res_wROI(:,3)=subset(vy,3);
res_wROI(:,4)=subset(vy,4);

figure
scatter(res_wROI(:,1), res_wROI(:,2),2), hold on

fprintf('\n -- Done. Selected ROI --\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:max(res_wROI(:,4))
     
    target=find(res_wROI(:,4)==i);
    
    if length(target)>min_length;
    
    tracks{count,1}(:,1)=res_wROI(target,3);    % frame
    tracks{count,1}(:,2)=res_wROI(target,1)*dx; % x in mum
    tracks{count,1}(:,3)=res_wROI(target,2)*dx; % y in mum
    
    count=count+1;
    else
    end
end

fprintf('\n -- Done. Generated tracks for @msdanalyzer --\n')



