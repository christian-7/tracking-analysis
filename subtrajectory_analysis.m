clear all, clc, close all
%% Load data and initialize parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt=0.5;                         % time step
dx=0.1;                         % pixel size ?m per pixel
segment=30;                     % Sm, segment length in frames
D=0.01;                        % D from first part of MSD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% new=dlmread('chaotic4.txt');
% 
% pos(:,1)=nonzeros(new(:,2)*dx);                                      % x coord in ?m
% pos(:,2)=nonzeros(new(:,3)*dx);                                      % y coord in ?m
% help=(new(:,1))+1;                                                   % helper variable
% pos(:,3)=help(1:length(pos),1);                                      % time in frames
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('tracks2.mat')

pos1(:,1)=nonzeros(tracks.Ch2(:,2)*dx);                                      % x coord in ?m
pos1(:,2)=nonzeros(tracks.Ch2(:,3)*dx);                                      % y coord in ?m
help=(tracks.Ch2(:,1))+1;                                                   % helper variable
pos1(:,3)=help(1:length(pos1),1);                                             % time in frames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%4%%%%%%%%%%%%%%%%

num_steps=length(pos1);
figure('Position',[200 400 900 500],'name','Overview Figure: scatter, velocity, MSD, cum. displacement')
h=gcf;
set(h,'PaperOrientation','landscape');

%% Select region of interest (time)

tmin=140;
tmax=230;

pos1(:,3)=pos1(:,3)-min(pos1(:,3))+1;

target=find(pos1(:,3)<tmax & pos1(:,3)>tmin);
pos(:,1)=pos1(target,1);
pos(:,2)=pos1(target,2);
pos(:,3)=pos1(target,3)-tmin;

num_steps=length(pos);


%% Plot Data

subplot(2,3,1)
line(pos(:,1),pos(:,2));hold on;
scatter(pos(:,1),pos(:,2),3,pos(:,3));hold on;
plot(pos(1,1),pos(1,2),'*b','MarkerSize',12);hold on;
text(pos(1,1),pos(1,2), 'Start');
plot(pos(length(pos),1),pos(length(pos),2),'+b','MarkerSize',12);hold on;
text(pos(length(pos),1),pos(length(pos),2),'End');hold on;
title('XY scatter trajectory');
xlabel('x (\mu m)','FontSize',12);
ylabel('y (\mu m)','FontSize',12);
% colorbar('northoutside');

%% Plot displacement from origin

dcum = zeros(num_steps,1);
dcum(1,1)=0;

for k=2:num_steps;
    
    d(k,1)=sqrt(((pos(k,1)-pos(1,1))^2)+((pos(k,2)-pos(1,2))^2));

    dcum(k)=d(k)+dcum(k-1);
end

subplot(2,3,4)
scatter(pos(:,3),d,2,pos(:,3));
title('Distance from origin');
xlabel('time (s)','FontSize',12);
ylabel('distance from origin (\mu m)','FontSize',12);
box on;

subplot(2,3,5)
scatter(pos(:,3),dcum,2,pos(:,3));
set(gca,'xscale','log')
set(gca,'yscale','log')
% axis([1 100 0.001 1e5])
title('Cumulative Distance from Origin');
xlabel('time (s)','FontSize',12);
ylabel('cumulative distance from origin (\mu m)','FontSize',12);
box on;

%% Calculate MSD (from Nathanael)

tab(:,1)=pos(:,3); %1:1:length(pos);
tab(:,2)=pos(:,1);
tab(:,3)=pos(:,2);
X=tab(:,2:3);

% dx=0.1;     % pixel size
% dt=0.5;       % time step

% X=data(:,2:3);

frame=tab(:,1);
frame=frame-min(frame);           % if it does not start with 0
% frame=frame/dt;                 % frame in sec
N=(max(frame)-min(frame)+1);      % number of frames

msd=zeros(1,N);
Time = [0:N-1]*dt;
f=zeros(1,N);


% Find the frames that have been recorded

for j=0:N
    if ~isempty(find(frame==j))
        f(j+1)=find(frame==j,1);
    end
end



for i=1:N-1
    c=0;
    
    for j=0:N-i-1
        if f(i+j+1)>0 && f(j+1)>0
            c=c+1;
            msd(i+1)=msd(i+1)+ sum((X(f(i+j+1),:)-X(f(j+1),:)).^2 );
        end
    end
    if c>0
        msd(i+1)=msd(i+1)/c;
    end
end


msd=msd/dt;



subplot(2,3,2)
box on;
plot(Time,msd,'-r');hold on;
title('Mean square displacement');
xlabel('time step','FontSize',12);
ylabel('MSD (\mu m^2)','FontSize',12);
% axis([0 50 0 5]);
% legend('Nathanael');
% set(leg2,'FontSize',12);
hold on;

% clear dt;

%% My MSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%transform time in frames

frame=[];
frame=pos(:,3); 
% frame=frame/dt;
% frame=frame-min(frame);
% frame(1,1)=1;

% i = frame --> Reihe
% j = gap; --> Spalte

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate MSD

msd=zeros(max(frame), 50);
msd2=zeros(50, 3);

for i=1:max(frame);    % find frame, for all frames
    vx=find(frame == i);
    
    if isempty(vx)==1; % if frame does not exist, skip   
    else
        
    
    for j=1:floor(max(frame)/4); % time gap
    
    if vx+j>length(pos) || i+j~=frame(vx+j) % if point plus gap exeeds lentgh, skip
    
        msd(i,j)=0;
    
    else
        
          msd(i,j)=((pos(vx,1)-pos(vx+j,1))^2)+((pos(vx,2)-pos(vx+j,2))^2);
    
%         msd2(j,1)=j;
%         msd2(j,2)=mean(nonzeros(msd(:,j)));
%         msd2(j,3)=std(nonzeros(msd(:,j)));
    end
    
    end
        
    end
      
    
end

for m=1:floor(max(frame)/4);
    
        msd2(m,1)=m;
        msd2(m,2)=mean(nonzeros(msd(:,m)));
        msd2(m,3)=std(nonzeros(msd(:,m)));
        
end

msd2(:,4)=msd2(:,2)/dt;
msd2(:,5)=msd2(:,1)*dt;
msd2(:,6)=msd2(:,3)/dt;

subplot(2,3,2)
box on;
% errorbar(msd2(:,5),msd2(:,4),msd2(:,6)); hold on;
plot(msd2(:,5),msd2(:,4),'-g'); hold on;
% plot(msd2(:,5), 4*D*msd2(:,5),'-g')



%% Velocity

vel=zeros(length(pos),1);

for k=2:num_steps;
    
    vel(k,1)=sqrt(((pos(k,1)-pos(k-1,1))^2)+((pos(k,2)-pos(k-1,2))^2));

    
end

vel=vel/dt; % vel in ?m/s
subplot(2,3,3)
scatter(pos(:,3),vel,3,pos(:,3));
title('Velocity');
xlabel('time (s)','FontSize',12);
ylabel('velocity (\mu m/s)','FontSize',12);
box on;


subplot(2,3,6)
hist(vel,20)
title('Histogram of Velocity');
xlabel('velocity (\mu m/s)','FontSize',12);
ylabel('counts','FontSize',12);
box on;

figure('Position',[200 400 300 300],'name','MSD for fitting')

scatter(msd2(2:15,5),msd2(2:15,4),'*r')
title('Mean square displacement');
xlabel('time step (s)','FontSize',12);
ylabel('MSD (\mu m^2)','FontSize',12);
box on

%% Calculate confinement

% generate variable frame

frame=pos(:,3);             % time step in seconds
frame=frame/dt;             % time step in frames
frame=frame-min(frame);     % starting from 0
frame(1,1)=1;               % starting from 1

% i = frame --> Reihe
% j = gap; --> Spalte

prob=[];    %zeros(max(frame), 5);
prob2=[];   %zeros(5, 2);
d=[];
vx=[];
vy=[];

c=1;

for i=1:max(frame)-segment;                         % for all frames
    vx=find(frame == i);                            % find frame i
    
    if isempty(vx)==1;                              % if frame does not exist, skip   
    else
        
     
    for j=4:segment;                                            % segment length
          
        vy=find(frame <= (i+j) & frame >= i );                  % select segment
        subset(:,1)=pos(vy);                                    % define segment as subset
        subset(:,2)=pos(vy,2);
        
        if length(vy)==1;                                       % if subset is only 1 frame --> distance is 0
                     R=0;
        else    
        
            for  k=2:length(subset);
                 d(k,1)=sqrt(((subset(k,1)-subset(1,1))^2)+((subset(k,2)-subset(1,2))^2));    % calculate the distance to each point in subset from point i  
            end
        R=max(d);                                                      % maximum distance within subset
        prob(i:(i+j),c)=0.2048-2.5117*((D*j)./(R^2));                  % probability within subset
%       prob(c,i:(i+j))=((D*j)./(R^2));

%       prob(c,i:(i+j))=horzcat(prob(c,i:(i+j)),((D*j)./(R^2)));
        c=c+1;  
        clear subset
        end
    
%     c=c+1;    
        
    end
    clear vx vy R d;
    
    end
   
end
clear subset

for l=1:length(frame)
    
prob2(l,1)=l;                           % frame
prob2(l,2)=mean(nonzeros(prob(l,:)));   % this is psi

end


L=[];

for i=1:length(prob2)
    
    if 10.^(prob2(i,2))>0.1
       L(i,1)=0;
       
    else
        
        L(i,1)=((prob2(i,2))*(-1)-1);
        
    end
L(i,2)=i;
    
end

integral=sum(L(:,1)); % calculate integral

figure('Position',[200 20 900 300], 'name','Probability psi, log(psi), confinement index L')

subplot(1,3,1)
plot(prob2(:,1)*dt,prob2(:,2));hold on;
xlabel('time (s)','FontSize',12);
ylabel('mean log(\psi)','FontSize',12);

subplot(1,3,2)
plot(prob2(:,1)*dt,10.^(prob2(:,2)));hold on;
xlabel('time (s)','FontSize',12);
ylabel('mean \psi','FontSize',12);

subplot(1,3,3)
plot(L(:,2)*dt, L(:,1)); hold on;      
xlabel('time (s)','FontSize',12);
ylabel('confinement index L','FontSize',12);

%% Figures to save

figure('Position',[200 800 300 300], 'name','confinement index L')
plot(L(:,2)*dt, L(:,1)); hold on;      
xlabel('time (s)','FontSize',12);
ylabel('confinement index L','FontSize',12);


figure('Position',[500 800 600 300], 'name','XY scatter and confinement index L')
h=gcf;
set(h,'PaperOrientation','landscape');

subplot(1,2,1)
line(pos(:,1),pos(:,2));hold on;
scatter(pos(:,1),pos(:,2),3,pos(:,3));hold on;
plot(pos(1,1),pos(1,2),'*b','MarkerSize',12);hold on;
text(pos(1,1),pos(1,2), 'Start');
plot(pos(length(pos),1),pos(length(pos),2),'+b','MarkerSize',12);hold on;
text(pos(length(pos),1),pos(length(pos),2),'End');hold on;
title('XY scatter');
xlabel('x (\mu m)','FontSize',12);
ylabel('y (\mu m)','FontSize',12);
box on;
% colorbar('northoutside');


subplot(1,2,2)
plot(L(:,2)*dt, L(:,1)); hold on;
scatter(L(:,2)*dt, L(:,1),15,pos(:,3),'filled')
% plot(L(:,2)*dt, s3,'-r');       
xlabel('time (s)','FontSize',12);
ylabel('confinement index L','FontSize',12);
title('Confinement index L');


figure('Position',[10 200 300 300])
h=gcf;
set(h,'PaperOrientation','landscape');
hist(integral,10);
title('Hist of integral')
xlabel('integral','FontSize',12);
ylabel('counts','FontSize',12);
