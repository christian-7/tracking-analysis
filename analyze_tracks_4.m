clear, clc, close all
%% Load data and initialize parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt=0.5;                         % time step
dx=0.1;                         % pixel size ?m per pixel
segment=30;                     % Sm, segment length in frames
D=0.01;                         % D from first part of MSD
threshold=10;                    % threshold for the dwell time analysis

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

pos(:,1)=nonzeros(tracks.Ch2(:,2)*dx);                                      % x coord in ?m
pos(:,2)=nonzeros(tracks.Ch2(:,3)*dx);                                      % y coord in ?m
help=(tracks.Ch2(:,1))+1;                                                   % helper variable
pos(:,3)=help(1:length(pos),1);                                             % time in frames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_steps=length(pos);
figure('Position',[200 400 900 500],'name','Overview Figure: scatter, velocity, MSD, cum. displacement')
h=gcf;
set(h,'PaperOrientation','landscape');

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
    
    dist(k,1)=sqrt(((pos(k,1)-pos(1,1))^2)+((pos(k,2)-pos(1,2))^2));

    dcum(k)=dist(k)+dcum(k-1);
end

subplot(2,3,4)
scatter(pos(:,3),dist,2,pos(:,3));
title('Distance from origin');
xlabel('time (s)','FontSize',12);
ylabel('distance from origin (\mu m)','FontSize',12);

subplot(2,3,5)
scatter(pos(:,3),dcum,2,pos(:,3));
set(gca,'xscale','log')
set(gca,'yscale','log')
% axis([1 100 0.001 1e5])
title('Cumulative Distance from Origin');
xlabel('time (s)','FontSize',12);
ylabel('cumulative distance from origin (\mu m)','FontSize',12);

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
%frame=frame/dt;                 % frame in sec
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


% msd=msd*dx^2/dt;
msd=msd/dt;


subplot(2,3,2)
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

msd2(:,4)=msd2(:,2)/dt; %dx^2/dt;           % mean in ?m2/sec
msd2(:,5)=msd2(:,1)*dt;                     % time in sec
msd2(:,6)=msd2(:,3)/dt; %/dx^2/dt;          % std in ?m2/s

subplot(2,3,2)
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


subplot(2,3,6)
hist(vel,20)
title('Histogram of Velocity');
xlabel('velocity (\mu m/s)','FontSize',12);
ylabel('counts','FontSize',12);

figure('Position',[200 400 300 300],'name','MSD for fitting')

scatter(msd2(1:6,5),msd2(1:6,4),'*r')
title('Mean square displacement');
xlabel('time step (s)','FontSize',12);
ylabel('MSD (\mu m^2)','FontSize',12);
box on

%% Calculate confinement

% generate variable frame

frame=pos(:,3);             % time step in seconds
frame=frame;%/dt;             % time step in frames
frame=frame-min(frame);     % starting from 0
frame=frame+1;
% frame(1,1)=1;               % starting from 1

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


figure('Position',[500 800 800 300], 'name','XY scatter and confinement index L')
h=gcf;
set(h,'PaperOrientation','landscape');

subplot(1,2,1)
line(pos(:,1),pos(:,2));hold on;
scatter(pos(:,1),pos(:,2),5,pos(:,3));hold on;
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
line([1 max(L(:,2)*dt)],[5 5],'Color','red')
axis([0 max(L(:,2)*dt) 0 70])
% plot(L(:,2)*dt, s3,'-r');       
xlabel('time (s)','FontSize',12);
ylabel('confinement index L','FontSize',12);
title('Confinement index L');

%% 

figure('Position',[500 800 500 600], 'name','Compare time with displacement and confinement')
h=gcf;
set(h,'PaperOrientation','portrait');

subplot(4,1,1)
line(pos(:,3)*dt,pos(:,1));hold on;
scatter(pos(:,3)*dt,pos(:,1),5,pos(:,3)*dt);hold on;
xlabel('time (s)','FontSize',12);
ylabel('x (\mu m)','FontSize',12);
box on;



subplot(4,1,2)
line(frame*dt,pos(:,2));hold on;
scatter(frame*dt,pos(:,2),5,pos(:,3));hold on;
xlabel('time (s)','FontSize',12);
ylabel('y (\mu m)','FontSize',12);
box on;

subplot(4,1,3)
plot(L(:,2)*dt, L(:,1)); hold on;
scatter(L(:,2)*dt, L(:,1),15,pos(:,3),'filled')     
xlabel('time (s)','FontSize',12);
ylabel('confinement index L','FontSize',12);
box on;

subplot(4,1,4)
plot(frame*dt,dist); hold on;
scatter(frame*dt,dist,2,pos(:,3)*dt);
xlabel('time (s)','FontSize',12);
ylabel('distance from origin (\mu m)','FontSize',12);
box on;

%% Identify Wells and measure dwell time

dwell=find(L(:,1) > threshold); % identify well with L>threshold

% dwell=dwell-min(dwell)+1; 

MSDwell=struct('well', [],'dwell_time', [],'coord_well',[],'MSD', [],'fit', []);

wells=[];
dwell_time=[];
p=1;

for o=1:length(dwell)-1;
    
    if dwell(o+1,1)==dwell(o,1)+1;
    
    wells(o,p)=dwell(o,1);               % frame that fall into the well column = well / row = frame
    MSDwell.well(o,p)=dwell(o,1);
    
    else
        
    p=p+1;
    wells(o,p)=dwell(o,1);
    MSDwell.well(o,p)=dwell(o,1);
    
    end
    
end

[m,n]=size(wells);

for q=1:n;
dwell_time(q,1)=length(nonzeros(wells(:,q)))
MSDwell.dwell_time(q,1)=length(nonzeros(wells(:,q)));
end

figure('Position',[500 800 400 300], 'name',['Dwell time in wells with Confinement Index > ' num2str(threshold)])
h=gcf;
set(h,'PaperOrientation','portrait');
hist(dwell_time)
xlabel('dwell time (s)','FontSize',12);
ylabel('counts','FontSize',12);
box on;

%% MSD per well
figure('Position',[500 800 600 300], 'name','MSD per well analysis using linear regression')
for r=1:n;

% MSDwell.coord_well(:,1)=pos(nonzeros(MSDwell.well(:,r)),1);
% MSDwell.coord_well(:,2)=pos(nonzeros(MSDwell.well(:,r)),2);       
% MSDwell.coord_well(:,3)=pos(nonzeros(MSDwell.well(:,r)),3);

frame=[];
pos3=[];

% pos3(:,1)=MSDwell.coord_well(:,1);
% pos3(:,2)=MSDwell.coord_well(:,2);
% pos3(:,3)=MSDwell.coord_well(:,3);

pos3(:,1)=pos(nonzeros(MSDwell.well(:,r)),1);
pos3(:,2)=pos(nonzeros(MSDwell.well(:,r)),2);
pos3(:,3)=pos(nonzeros(MSDwell.well(:,r)),3);



frame=pos3(:,3); 
frame=frame-min(frame)+1;
% frame=frame*dt;


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
    
    if vx+j>length(pos3) || i+j~=frame(vx+j)             % if point plus gap exeeds length, skip
    
        msd(i,j)=0;
    
    else
        
          msd(i,j)=((pos3(vx,1)-pos3(vx+j,1))^2)+((pos3(vx,2)-pos3(vx+j,2))^2);
    
%         msd2(j,1)=j;
%         msd2(j,2)=mean(nonzeros(msd(:,j)));
%         msd2(j,3)=std(nonzeros(msd(:,j)));
    end
    
    end
        
    end
      
    
end

        MSDwell.MSD{r,1}(1,1)=0;
        MSDwell.MSD{r,1}(1,2)=0;
        MSDwell.MSD{r,1}(1,3)=0;

for m=1:floor(max(frame)/4);
    
        MSDwell.MSD{r,1}(m+1,1)=m;
        MSDwell.MSD{r,1}(m+1,2)=mean(nonzeros(msd(:,m)));
        MSDwell.MSD{r,1}(m+1,3)=std(nonzeros(msd(:,m)));
        
end

MSDwell.MSD{r,1}(:,4)=MSDwell.MSD{r,1}(:,2)/dt; % dx^2/dt;           % mean in ?m2/sec
MSDwell.MSD{r,1}(:,5)=MSDwell.MSD{r,1}(:,1)*dt;                      % time in sec
MSDwell.MSD{r,1}(:,6)=MSDwell.MSD{r,1}(:,3)/dt; % dx^2/dt;           % std in ?m2/s

% Perfom linear fitting

a1 = polyfit(MSDwell.MSD{r,1}(1:5,5),MSDwell.MSD{r,1}(1:5,4),1);

MSDwell.fit{r,1}(:,1)=a1(:,1);
MSDwell.fit{r,1}(:,2)=a1(:,2);
xfit=0:0.1:2.5;

subplot(1,2,1);
scatter(MSDwell.MSD{r,1}(1:5,5),MSDwell.MSD{r,1}(1:5,4),'*'); hold on;
plot(xfit,a1(1)*xfit+a1(2));hold on;
legendInfo{r} = ['MSD well = ' num2str(r)]; 


% scatter(MSDwell.MSD{r,1}(1:floor(max(frame)/4),5),MSDwell.MSD{r,1}(1:floor(max(frame)/4),4),'*'); hold on;
% a1 = polyfit(MSDwell.MSD{r,1}(1:floor(max(frame)/4),5),MSDwell.MSD{r,1}(1:floor(max(frame)/4),4),1);


% a2=fit(x,y,'r_c^2*(1-A1*exp((-4*A2*A3*x)/r_c^2))','StartPoint', [0.01, 2, 5])

subplot(1,2,2);
scatter(MSDwell.dwell_time(r,1),MSDwell.fit{r,1}(1,1));hold on;
xlabel('dwell time (s)','FontSize',12);
ylabel('D_{conf} (\mum^2/s)','FontSize',12);

clear a1;

end
legend(legendInfo);



%% Perform least square fitting to determine D_conf

% a2=fit(x,y,'r_c^2*(1-A1*exp((-4*A2*A3*x)/r_c^2))','StartPoint', [0.1,0.1,0.001,0.01] )
% figure
% plot(a2,x,y);hold on;
% scatter(x,y);
figure('Position',[500 800 600 600], 'name','MSD per well analysis using conf model')

for r=1:n;

    x=MSDwell.MSD{r,1}(1:6,5);
    y=MSDwell.MSD{r,1}(1:6,4);
    
xfit=0:0.1:2.5;

fun = @(p,x) p(1)^2*(1-p(2)*exp((-4*p(3)*p(4)*x)/p(1)^2));
% fun = @(p,x) p(1)^2*(1-1*exp((-4*2*p(2)*x)/p(1)^2));


pguess = [0.05,1,2,0.001]; % r, A1, A2, D
% pguess = [0.05,0.001]; % r, D
[p,fminres] = lsqcurvefit(fun,pguess,x,y)

MSDwell.fit{1,1}(r,1)=p(1); % radius
MSDwell.fit{1,1}(r,2)=p(2); % A1, D
MSDwell.fit{1,1}(r,3)=p(3); % A2
MSDwell.fit{1,1}(r,4)=p(4); % D
MSDwell.fit{1,1}(r,5)=fminres; % minres

subplot(2,2,1)
plot(xfit,p(1)^2*(1-p(2)*exp((-4*p(3)*p(4)*xfit)/p(1)^2)));hold on;
% plot(xfit,p(1)^2*(1-1*exp((-4*2*p(2)*xfit)/p(1)^2)));hold on;
scatter(x,y);
xlabel('lag time (s)','FontSize',12);
ylabel('MSD (\mum^2/s)','FontSize',12);

clear x y;
end

subplot(2,2,2)
scatter(MSDwell.dwell_time,MSDwell.fit{1,1}(:,4)) %4
xlabel('dwell time (s)','FontSize',12);
ylabel('D_{conf} (\mum^2/s)','FontSize',12);

subplot(2,2,3)
scatter(MSDwell.fit{1,1}(:,1),MSDwell.fit{1,1}(:,4)) %4
xlabel('radius (\mum)','FontSize',12);
ylabel('D_{conf} (\mum^2/s)','FontSize',12);

subplot(2,2,4)
scatter(MSDwell.dwell_time,MSDwell.fit{1,1}(:,1))
xlabel('dwell time (s)','FontSize',12);
ylabel('radius (\mum)','FontSize',12);
