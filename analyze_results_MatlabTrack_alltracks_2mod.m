clear all, close all, clear all, clc

%% Read file --> Output from readUtrackTracks.py

filename='tracksFinal_tracksXYT.txt';

delimiterIn = '\t';
headerlinesIn = 1;
tracks = importdata(filename,delimiterIn,headerlinesIn);

col1=tracks.textdata(2:end,1);

for index=1:length(col1);
    
    minus=find(col1{index,1}=='-');
    cell2mat(col1(index,1));
    trackID(index,1)=str2num(ans(1:(minus-1)));
    index=index+1;
end

tracks=tracks.data;
tracks(:,8)=trackID; % put the track ID in column 8

clear index ans col1 delimiterIn headerlinesIn trackID minus;

%% Put all individual trajectories in structure indTracks
%  create new structure inTracks

indTracks=struct('tracks',[], 'gaps', [], 'gapframes', [], 'offtime', [],'blinkframes', [], 'ontime' ,[], 'nbrofgaps' ,[], 'gaplength' ,[],'nbrofblinks' ,[], 'blinklength' ,[], 'pdist',[]);

for number=1:max(tracks(:,8));

    target=find(tracks(:,8)==number);
    
    indTracks.tracks{number,1}=tracks(target);          % Track ID, i.e. number of trajectory
    indTracks.tracks{number,1}(:,2)=tracks(target,2);   % Frame
    indTracks.tracks{number,1}(:,3)=tracks(target,3);   % x position
    indTracks.tracks{number,1}(:,4)=tracks(target,4);   % Std X
    indTracks.tracks{number,1}(:,5)=tracks(target,5);   % y position
    indTracks.tracks{number,1}(:,6)=tracks(target,6);   % Std Y
    indTracks.tracks{number,1}(:,7)=tracks(target,7);   % Amplitude
        
end

 
%%  Find gaps (set1) and blinks (set0) in X coordinate

for track=1:length(indTracks.tracks);

vx=isnan(indTracks.tracks{track,1}(:,2));  % finds nan in subset x coordinate
                                           % nan (gaps) --> 1
for index=1:length(vx); 
    
        if  vx(index)==1;                                                           % if nan, copy  frame number in new
            indTracks.gaps{track,1}(index,1)=indTracks.tracks{track,1}(index,1);    % new(:,1) --> gap frames
           
        else indTracks.gaps{track,1}(index,2)=indTracks.tracks{track,1}(index,1);   % new(:,2) --> blink frames
       
        end
end



indTracks.gapframes(track,1)=length(nonzeros(indTracks.gaps{track,1}(:,1)));                             % number of gap frames
indTracks.offtime(track,1)=indTracks.gapframes(track,1)/length(indTracks.gaps{track,1}(:,1));            % percentage off time

if  max(indTracks.gaps{track,1}(:,1))==0;
    indTracks.blinkframes(track,1)=1;
else
indTracks.blinkframes(track,1)=length(nonzeros(indTracks.gaps{track,1}(:,2)));                           % number of blink frames
end

indTracks.ontime(track,1)=indTracks.blinkframes(track,1)/length(indTracks.gaps{track,1}(:,1));           % fraction of time the molecue stays on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear track;

%% Analyze Gaps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for track=1:length(indTracks.gaps);

blinks=nonzeros(indTracks.gaps{track,1}(:,1)); % isolate gap frames 

blinks2=[];
blinks3=[];
gaps2=[];

% find gap start

if          isempty(blinks) == 1;               % if there is no gap > blinks is empty > set nbr to 0
            indTracks.nbrofgaps(track,1)=0;
            %indTracks.gaplength(track,1)=0;
else            

for index=1:length(blinks)-1; 
    
        if  blinks(index+1) == blinks(index)+1;
            %gaps(index,2)=gaps(index,1);

            
        else blinks2(index,1)= blinks(index,1)+1;
             %gaps(index,2)=gaps(index,1);
        end
end

% find gap end

for index=2:length(blinks); 
    
        if  blinks(index-1) == blinks(index)-1;
            %gaps(index,2)=gaps(index,1);
             
        else blinks3(index,1)= blinks(index,1)-1;
             %gaps(index,2)=gaps(index,1);
        end
end

if  isempty(blinks3) == 1;
    indTracks.nbrofgaps(track,1)=1;   % number of gaps
    
else 
    
gaps2(:,1)=nonzeros(blinks3(:,1));                               % gap start
gaps2(:,2)=nonzeros(blinks2(:,1));                               % gap stop
indTracks.nbrofgaps(track,1)=length(nonzeros(gaps2(:,1)))+0;     % number of gaps

end

% calculate gap length

if          indTracks.nbrofgaps(track,1)==1;                    % if the number of gaps is 1 > the length is equal to the length of blinks
            indTracks.gaplength{track,1}=length(blinks);
 
else
            
if isempty(gaps2)==1;                                           % if there is no gap start > there is no gap > gaplength = 1
   indTracks.gaplength{track,1}=0;
   
   
else    for     index=1:length(gaps2(:,1));
                indTracks.gaplength{track,1}(index,1)=gaps2(index,1)-gaps2(index,2)+1;  % stop - start frame  
        end
end

clear blinks gaps2 blinks2 blinks3
 
end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyze Blinks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% isolate gap frames to calculate blinks

for track=1:length(indTracks.gaps);

blinks=nonzeros(indTracks.gaps{track,1}(:,2)); % isolate blink frames to calculate gaps

blinks2=[];
blinks3=[];
gaps2=[];

% find blink start

for index=1:length(blinks)-1; 
    
        if  blinks(index+1) == blinks(index)+1;
            %gaps(index,2)=gaps(index,1);
             
        else blinks2(index,1)=blinks(index,1)+1;
             %gaps(index,2)=gaps(index,1);
        end
end

% find blink end

for index=2:length(blinks); 
    
        if  blinks(index-1) == blinks(index)-1;
            %gaps(index,2)=gaps(index,1);
             
        else blinks3(index,1)=blinks(index,1)-1;
             %gaps(index,2)=gaps(index,1);
        end
end

if  isempty(blinks3) == 1;
    indTracks.nbrofblinks(track,1)=1;                              % number of gaps

else 
    
gaps2(:,1)=nonzeros(blinks3(:,1));                               % blink start
gaps2(:,2)=nonzeros(blinks2(:,1));                               % blink stop
indTracks.nbrofblinks(track,1)=length(nonzeros(gaps2(:,1)))+1;   % number of blinks

end

% calculate blink length

if isempty(gaps2)==1;
   indTracks.blinklength{track,1}=0;
   
else    for     index=1:length(gaps2(:,1));
                indTracks.blinklength{track,1}(index,1)=gaps2(index,1)-gaps2(index,2)+1;  % stop - start frame  
        end
end

clear blinks gaps2 blinks2 blinks3
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Data Visualization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Show histograms of track length, blink length, gap length

figure('Position',[100 500 700 200])

% Track length

for index=1:max(tracks(:,8));
    
    target=find(tracks(:,8) == index);
    track_target=tracks(target);
    tracklength(:,index)=length(track_target);
    
end

subplot(1,3,1)
hist(tracklength,20);  
title('Total track length (frames)')
xlabel('track length (frames)');
ylabel('counts');
  
% Blink length
  
  blinklengthHist=[];

for index=1:length(indTracks.blinklength)
    
    blinklengthHist=vertcat(blinklengthHist, nonzeros(indTracks.blinklength{index,1}));
end

subplot(1,3,2)
hist(blinklengthHist,20)
title('blink length')
xlabel('blink time (frames)');
ylabel('counts');


% Gap Length

    gaplengthHist=[];

for index=1:length(indTracks.gaplength)
    
    gaplengthHist=vertcat(gaplengthHist, nonzeros(indTracks.gaplength{index,1}));
end

subplot(1,3,3)
hist(gaplengthHist,20)
title('Gap length')
xlabel('dark time (frames)');
ylabel('counts');

% Plot histograms from frame analysis 

figure('Position',[100 100 700 500])

subplot(2,3,1)
hist(indTracks.gapframes,20)
title('# Gap Frames')
xlabel('number of gap frames');
ylabel('counts');

subplot(2,3,4)
hist(indTracks.nbrofgaps,20)
title('# of gaps')
xlabel('number of gaps');
ylabel('counts');

subplot(2,3,2)
hist(indTracks.blinkframes,20)
title('# Blink Frames')
xlabel('number of blink frames');
ylabel('counts');

subplot(2,3,5)
hist(indTracks.nbrofblinks,20)
title('# of blinks')
xlabel('number of blinks');
ylabel('counts');

subplot(2,3,3)
hist(indTracks.offtime,20)
title('fraction off time')
xlabel('fraction off time');
ylabel('counts');

subplot(2,3,6)
hist(indTracks.ontime,20)       % fraction of time the molecue stays on
title('fraction on time')
xlabel('fration on time');
ylabel('counts');


%% Normalized integrated counts for gap length and number of blinks

binrange=0:max(gaplengthHist);
[bincounts] = histc(gaplengthHist,binrange);

figure('Position',[100 500 700 500])

subplot(2,2,1)
bar(binrange,bincounts,'histc')
title('dark time (frames)')
xlabel('dark time (frames)');
ylabel('counts');

for index=2:(length(bincounts));
    comcounts(1,1)=bincounts(1);
    comcounts(index)=bincounts(index)+comcounts(index-1);    
end

comcountsNorm=comcounts/max(comcounts);

subplot(2,2,2)
plot(binrange, comcountsNorm);
title('dark time (frames)');
xlabel('dark time (frames)');
ylabel('normalized integrated counts');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

binrangeNbr=0:max(indTracks.nbrofblinks);
[bincountsNbr] = histc(indTracks.nbrofblinks,binrangeNbr);

subplot(2,2,3)
bar(binrangeNbr,bincountsNbr,'histc')
title('number of blinks')
xlabel('number of blinks');
ylabel('counts');

for index=2:(length(bincountsNbr));
    comcountsNbr(1,1)=bincountsNbr(1);
    comcountsNbr(index)=bincountsNbr(index)+comcountsNbr(index-1);    
end

comcountsNbrNorm=comcountsNbr/max(comcountsNbr);

subplot(2,2,4)
plot(binrangeNbr, comcountsNbrNorm);
title('number of blinks');
xlabel('number of blinks');
ylabel('normalized integrated counts');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Individual tracks%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% find xy coordinates for all tracks and copy them in indTracks.pdist
%  Calculate spatial scattering

for track=1:length(indTracks.tracks);

vx=isnan(indTracks.tracks{track,1}(:,2));      % finds nan in subset x coordinate
                                                % nan --> 1
                                                
if length(vx)>2                                                
                                                
clusterx=[];
clustery=[];
                                               
for index=1:length(vx); 
    
   
        
        if  vx(index)==0;                                                             % if nan, copy  frame number in new

            clusterx=vertcat(clusterx,indTracks.tracks{track,1}(index,2));
            clustery=vertcat(clustery,indTracks.tracks{track,1}(index,4));
 
        else  1;  
       
        end
end

            cluster(:,1)=clusterx;
            cluster(:,2)=clustery;
            
            indTracks.pdist{track,1}=cluster;
            indTracks.pdist{track,2}=pdist(cluster,'euclidean');        % distance distribution
            indTracks.pdist{track,3}=mean(indTracks.pdist{track,2});    % mean distance
            indTracks.pdist{track,4}=max(indTracks.pdist{track,2});     % max distance
            
clear cluster clusterx clustery            

else        indTracks.pdist{track,1}=0;
            indTracks.pdist{track,2}=0;        % distance distribution
            indTracks.pdist{track,3}=0;        % mean distance
            indTracks.pdist{track,4}=0;        % max distance
        
 end

end

figure('Position',[100 500 700 300])

subplot(1,2,1)
hist(nonzeros(cell2mat(indTracks.pdist(:,3))),20);
title('Hist of mean distance');
xlabel('distance (pxl)');
ylabel('counts');

subplot(1,2,2)
hist(nonzeros(cell2mat(indTracks.pdist(:,4))), 20);
title('Hist of max distance');
xlabel('distance (pxl)');
ylabel('counts');

%% PLot individual molecule scatters in subplots

max=50;
figure('Position',[100 500 700 300])

m=floor(sqrt(max))+1;

for track=1:max;

vx=isnan(indTracks.tracks{track,1}(:,2));      % finds nan in subset x coordinate
                                                % nan --> 1
                                                
%if length(vx)>2                                                
                                                
clusterx=[];
clustery=[];
                                               
for index=1:length(vx); 
    
   
        
        if  vx(index)==0;                                                             % if nan, copy  frame number in new

            clusterx=vertcat(clusterx,indTracks.tracks{track,1}(index,2));
            clustery=vertcat(clustery,indTracks.tracks{track,1}(index,4));
 
        else  1;  
       
        end
end


subplot(m, m, track); hold on
scatter(clusterx, clustery)

%clear clusterx clustery
             

end

%% PLot individual molecule scatters in the same plot (1) and show histograms of xy dimension, mean and max

max=200;
figure('Position',[100 500 700 500])
m=floor(sqrt(max))+1;
allclustersx=[];
allclustersy=[];

for track=1:max;

vx=isnan(indTracks.tracks{track,1}(:,2));      % finds nan in subset x coordinate
                                                % nan --> 1
                                                
%if length(vx)>2                                                
                                                
clusterx=[];
clustery=[];
                                               
for index=1:length(vx); 
    
   
        
        if  vx(index)==0;                                                             % if nan, copy  frame number in new

            clusterx=vertcat(clusterx,indTracks.tracks{track,1}(index,2));
            clustery=vertcat(clustery,indTracks.tracks{track,1}(index,4));
 
        else  1;  
       
        end
end

clusterx=clusterx-min(clusterx);
clustery=clustery-min(clustery);

allclustersx=vertcat(allclustersx,clusterx);
allclustersy=vertcat(allclustersy,clustery);

subplot(2,3,1)
scatter(clusterx, clustery); hold on;
title('Scatter of all molecules');
xlabel('distance (pxl)');
ylabel('distance (pxl)');


clear clusterx clustery
% else             
% end
end

subplot(2,3,2)
hist(allclustersx,20)
title('X diameter');
xlabel('distance (pxl)');
ylabel('counts');

subplot(2,3,3)
hist(allclustersy,20)
title('Y diameter');
xlabel('distance (pxl)');
ylabel('counts');

subplot(2,3,5)
hist(nonzeros(cell2mat(indTracks.pdist(:,3))),20);
title('Mean point distance');
xlabel('distance (pxl)');
ylabel('counts');

subplot(2,3,6)
hist(nonzeros(cell2mat(indTracks.pdist(:,4))), 20);
title('Max point distance');
xlabel('distance (pxl)');
ylabel('counts');

%% Generate normalized histograms of xy width and mean/max distance

[f,x]=hist(allclustersx,20); 
subplot(2,2,1)
bar(x,f/sum(f));
title('Hist of x diameter');
xlabel('distance (pxl)');
ylabel('norm counts');

clear f x

[f,x]=hist(allclustersy,20); 
subplot(2,2,2)
bar(x,f/sum(f));
title('Hist of y diameter');
xlabel('distance (pxl)');
ylabel('norm counts');
clear f x

[f,x]=hist(nonzeros(cell2mat(indTracks.pdist(:,3))),20); 
subplot(2,2,3)
bar(x,f/sum(f));
title('Hist of mean distance');
xlabel('distance (pxl)');
ylabel('norm counts');
clear f x

[f,x]=hist(nonzeros(cell2mat(indTracks.pdist(:,4))),20); 
subplot(2,2,4)
bar(x,f/sum(f));
title('Hist of max distance');
xlabel('distance (pxl)');
ylabel('norm counts');

%% Calculate PCF for x and y dimension

figure('Position',[100 500 700 300])

nbins=50;
[f,x]=hist(nonzeros(allclustersx),nbins); 

N=length(nonzeros(allclustersx));
w=0.8038/nbins;

for index=1:length(f);

PCF1(index)=f(index)/(N*((x(index)+w).^2)-x(index).^2);

end

nbins=50;
[f,x]=hist(allclustersx,nbins); 

N=length(allclustersx);
w=0.8038/nbins;

for index=1:length(f);

PCF2(index)=f(index)/(N*((x(index)+w).^2)-x(index).^2);

end

subplot(1,2,1)
plot(x,PCF1,'blue',x,PCF2,'red');hold on;
title('PCF along x dimension');
xlabel('radius (pxl)');
ylabel('g(r)');
legend('nonzeros','all')

clear PCF1 PCF2 f x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbins=50;
[f,x]=hist(nonzeros(allclustersy),nbins); 

N=length(nonzeros(allclustersy));
w=0.8038/nbins;

for index=1:length(f);

PCF1(index)=f(index)/(N*((x(index)+w).^2)-x(index).^2);

end

clear f x

[f,x]=hist(allclustersy,nbins); 

N=length(allclustersy);
w=1.264/nbins;

for index=1:length(f);

PCF2(index)=f(index)/(N*((x(index)+w).^2)-x(index).^2);

end

subplot(1,2,2)
plot(x,PCF1,'blue',x,PCF2,'red');
title('PCF along y dimension');
xlabel('radius (pxl)');
ylabel('g(r)');
legend('nonzeros','all')

%% %% Normalize each cluster to its center of mass, plot 2D Histogram

max=length(indTracks.tracks);
figure('Position',[100 500 800 800])

m=floor(sqrt(max))+1;
allclustersCx=[];
allclustersCy=[];

for track=1:max;

vx=isnan(indTracks.tracks{track,1}(:,2));      % finds nan in subset x coordinate
                                                % nan --> 1
                                                
%if length(vx)>2                                                
                                                
clusterx=[];
clustery=[];
                                               
for index=1:length(vx); 
    
        if  vx(index)==0;                                                             % if nan, copy  frame number in new

            clusterx=vertcat(clusterx,indTracks.tracks{track,1}(index,2));
            clustery=vertcat(clustery,indTracks.tracks{track,1}(index,4));
 
        else  1;  
       
        end
end

clusterxC=sum(clusterx)/length(clusterx);
clusteryC=sum(clustery)/length(clustery);
clusterx=clusterx-clusterxC;
clustery=clustery-clusteryC;

allclustersCx=vertcat(allclustersCx,clusterx);
allclustersCy=vertcat(allclustersCy,clustery);

clear clusterx clustery
             
end

c=hist3([allclustersCx*100, allclustersCy*100],[20 20]);

subplot(2,2,1)
scatter(allclustersCx, allclustersCy)
subplot(2,2,2)
% surf(c)
hist3([allclustersCx*100, allclustersCy*100],[20 20])
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
%colorbar
subplot(2,2,3)
histfit(allclustersCx, 20)
subplot(2,2,4)
histfit(allclustersCy, 20)

%% Calculate gaussian PDF of x and y dimensions, overlay with histogram

figure('Position',[100 500 1000 300])

pdx=fitdist(allclustersCx(:),'normal')
pdy=fitdist(allclustersCy(:),'normal')

xPDF=-3:0.05:3;  
PDFx=pdf(pdx,xPDF); %PDF is a vector of y values: it's our fit
PDFy=pdf(pdy,xPDF);

% Now overlay the fit

 PDFx=PDFx/sum(PDFx); % normalise
 PDFy=PDFy/sum(PDFy);

% scale to y axis

y=ylim;
PDFx=PDFx*y(2);
PDFy=PDFy*y(2);

[f,x]=hist(allclustersCx,20); 
subplot(1,2,1)
bar(x*100,f/sum(f));
hold on
plot(xPDF*100,PDFx,'r-','LineWidth',2)
axis([-100 100 0 0.35])
title('PDF x dimension');
xlabel('radius (pxl)');
ylabel('pdf');
hold off

clear f x

[f,x]=hist(allclustersCy,20); 
subplot(1,2,2)
bar(x*100,f/sum(f));
hold on
plot(xPDF*100,PDFy,'r-','LineWidth',2)
axis([-100 100 0 0.35])
title('PDF y dimension');
xlabel('radius (pxl)');
ylabel('pdf');
hold on


%
%% 

clusterxC=sum(clusterx)/length(clusterx);
clusteryC=sum(clustery)/length(clustery);
clusterx=clusterx-clusterxC;
clustery=clustery-clusteryC;
scatter(clusterx, clustery)
clusterx=clusterx*100;
clustery=clustery*100;
scatter(clusterx, clustery)
axis([-40 40 -40 40])