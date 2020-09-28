change first line
clear;close all;clc;
eps=0.1;

Xsol= nan(4, 5001);
Xsol(1,1)=0.5;
Xsol(2, 1)=0.2;
Xsol(3, 1)=0.4;
Xsol(4, 1)=0.3;

tmax=500;
dt=0.1;
tspan=0:dt:tmax;
nT= length(tspan);
tspike1=[];
tspike2=[];

for i=1:nT-1
    X=Xsol(:, i);
    x1= X(1);
    x2= X(2);
    s1= X(3);
    s2= X(4);
    %euler's method
    dx1= (-x1*(1-x1)/10+0.5*s2+0.03)*dt;
    dx2= (-x2*(1-x2)/10+0.5*s1+0.03)*dt;
    ds1=-s1*0.1*dt;
    ds2=-s2*0.1*dt;
    x1new=x1+dx1;
    x2new=x2+dx2;
    s1new=s1+ds1;
    s2new=s2+ds2;
    if(abs(x1new-1)<0.1)
        x1new=0;
        s1new=s1new+eps;
        tspike1=[tspike1, (i-1)*dt];
    end
    if(abs(x2new-1)<0.1)
        x2new=0;
        s2new=s2new+eps;
        tspike2=[tspike2, (i-1)*dt];
    end
    Xsol(:,i+1)=[x1new, x2new, s1new, s2new];
end
%% new figure
fig = figure;
set(fig, 'position', get(0,'ScreenSize')); % Fullscreen
%% figure1
subplot(2,2,1)
plot(tspan, Xsol(1,:))
hold on
plot(tspan, Xsol(2,:))
y = Xsol(1:2,:); % y data
fSz = 20; % fontsize
myProcess(y,fSz); % postprocess

%% figure2
subplot(2,2,2)
plot(tspan, Xsol(3,:))
hold on
plot(tspan, Xsol(4,:))
y = Xsol(3:4,:); % y data
myProcess(y,fSz); % postprocess
%rasterplot(tspike1,10, 500,0)
%hold on
%rasterplot(tspike2, 10, 500,1)
%X=[x1;x2;s1;s2]

Xsoll= nan(4, 5001);
Xsoll(1,1)=0.5;
Xsoll(2, 1)=0.2;
Xsoll(3, 1)=0.4;
Xsoll(4, 1)=0.3;
tmax=1000;
dt=0.1;
tspan=0:dt:tmax;
nT= length(tspan);
tspike1=[];
tspike2=[];
for i=1:nT-1
    X=Xsoll(:, i);
    x1= X(1);
    x2= X(2);
    s1= X(3);
    s2= X(4);
    %euler's method
    dx1= (-x1*(1-x1)/10-0.1*s2+0.03)*dt;
    dx2= (-x2*(1-x2)/10-0.1*s1+0.03)*dt;
    ds1=-s1*0.1*dt;
    ds2=-s2*0.1*dt;
    x1new=x1+dx1;
    x2new=x2+dx2;
    s1new=s1+ds1;
    s2new=s2+ds2;
    
    if(abs(x1new-1)<0.1)
        x1new=0;
        s1new=s1new+eps;
        tspike1=[tspike1, (i-1)*dt];
    end
    if(abs(x2new-1)<0.1)
        x2new=0;
        s2new=s2new+eps;
        tspike2=[tspike2, (i-1)*dt];
    end
    Xsoll(:,i+1)=[x1new, x2new, s1new, s2new];
end
%% figure3
subplot(2,2,3);
plot(tspan, Xsoll(1,:))
hold on
plot(tspan, Xsoll(2,:))
y = Xsol(1:2,:); % y data
myProcess(y,fSz); % postprocess
%% figure4
subplot(2,2,4)
plot(tspan, Xsoll(3,:))
hold on
plot(tspan, Xsoll(4,:))
y = Xsol(3:4,:); % y data
myProcess(y,fSz); % postprocess
%rasterplot(tspike1,10, 500,0)
%hold on
%rasterplot(tspike2, 10, 500,1)
%% function
function myProcess(y,fSz)
y = reshape(y,1,[]);
set(gca,'FontSize',fSz);
ylabel('ylabel1','Interpreter','latex'); % ylabel
xlabel('xlabel1','Interpreter','latex'); % xlabel
set(gca,'box','on'); % box on
miny = min(y); % min y
maxy = max(y); % max y
ylim([miny maxy]); % ylim
set(gca, 'YTick', linspace(miny,maxy,3)); % YTick
set(gca,'YTickLabel',{'min','mid','max'}); % YTickLabel
grid on;
end

function rasterplot(times,numtrials,triallen, expectHeight)
%%%%%%%%%%%%%% Plot variables %%%%%%%%%%%%%%
plotwidth=1;     % spike thickness
plotcolor='k';   % spike color
trialgap=1.5;    % distance between trials
defaultfs=1000;  % default sampling rate
showtimescale=1; % display timescale
showlabels=1;    % display x and y labels
%%%%%%%%% Code Begins %%%%%%%%%%%%
%switch nin
% case 3 %no handle so plot in a separate figure
%figure
hresp=gca;
fs=defaultfs;
%case 4 %handle supplied
% hresp=varargin{1};
% if (~ishandle(hresp))
% error('Invalid handle');
% end
% fs=defaultfs;
% case 5 %fs supplied
% hresp=varargin{1};
%if (~ishandle(hresp))
%   error('Invalid handle');
% end
% fs = varargin{2};
% otherwise
% error ('Invalid Arguments');
%end
% plot spikes
trials=ceil(times/triallen);
reltimes=mod(times,triallen);
reltimes(~reltimes)=triallen;
numspikes=length(times);
xx=ones(3*numspikes,1)*nan;
yy=ones(3*numspikes,1)*nan;
yy(1:3:3*numspikes)=(trials-1)*trialgap+expectHeight;
yy(2:3:3*numspikes)=yy(1:3:3*numspikes)+1;
%scale the time axis to ms
xx(1:3:3*numspikes)=reltimes*1000/fs;
xx(2:3:3*numspikes)=reltimes*1000/fs;
xlim=[1,triallen*1000/fs];
axes(hresp);
h=plot(xx, yy, plotcolor, 'linewidth',2);
axis([xlim,0,(numtrials)*1.5]);
if (showtimescale)
    set(hresp, 'ytick', [],'tickdir','out');
else
    set(hresp,'ytick',[],'xtick',[]);
end
xlabel('Time(ms)');
ylabel('Trials');
end
