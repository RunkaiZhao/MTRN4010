function []=MTRN4010_FuzzyCtrl_()
% fuzzy control
% sim car movement, virtual tracking
rng('default');

clc; clear all; close all; dbstop if error;
set(0,'defaultaxesfontname','times new roman');

field.range=20;

% virtual target, you have to examine the effect of virtual target distance
Virtual=150; 

time.dt=2; time.T=1000;

[fig]=FigureNew(field);
[car]=CarNew(field); 
[car]=CarNow(car,time,0,0); 
[car]=CarShow(fig,car,0);
[target]=TargetNew(car);

% ou have to modify the fuzzy memberships for better results
fis_vel=readfis('MTRN4010_vel_4.fis');
fis_ang=readfis('MTRN4010_ang.fis');

for t=0:time.dt:time.T
  [target]=TargetNow(t,time,car,target,Virtual);
  [target]=TargetShow(fig,target,t);
  [ds]=FindDistance(car,target);
  vel=evalfis(ds,fis_vel);
  [dq]=FindAngular(car,target);
  ang=evalfis(dq,fis_ang);
  [car]=CarNow(car,time,vel,ang);
  [car]=CarShow(fig,car,t);
  disp(ds);
  
  if ds < 0.5
      break;
  end
  
end;

CarAxes(car);


function [fig]=FigureNew(field)
fig.fig=figure('units','normalized','position',[0.3 0.2 0.4 0.5]);
axis([-1 1 -1 1]*field.range*1.5); hold on; grid on; axis equal;
xlabel('x-direction'); ylabel('y-direction');
fig.ax=axis;

function [target]=TargetNew(car)
target.x=-car.x;
target.y=-car.y;
target.q=AngleWrap(rand*2*pi);
target.x_=target.x; target.y_=target.y; target.q_=target.q;
target.shape=[ 2 0; 1 1; -1 1; -1 -1; 1 -1; 2 0]';
target.hdL.shape=plot(target.shape(1,:),target.shape(2,:),...
  'color','r','linewidth',2);
target.hdL.shape_=plot(target.shape(1,:),target.shape(2,:),...
  'color','m','linewidth',2);
Rz=[  cos(target.q) -sin(target.q); 
      sin(target.q)  cos(target.q)];
shape=Rz*target.shape+repmat([target.x;target.y],1,6);
set(target.hdL.shape,'xdata',shape(1,:),'ydata',shape(2,:)); 
set(target.hdL.shape_,'xdata',shape(1,:),'ydata',shape(2,:)); 

function [target]=TargetNow(t,time,car,target,Virtual)
ds=Virtual*(1-t/time.T)^3;
target.x=target.x_-ds*cos(target.q);
target.y=target.y_-ds*sin(target.q);

function [target]=TargetShow(fig,target,t)
Rz=[  cos(target.q) -sin(target.q); 
      sin(target.q)  cos(target.q)];
shape=Rz*target.shape+repmat([target.x;target.y],1,6);
set(target.hdL.shape_,'xdata',shape(1,:),'ydata',shape(2,:)); 
title(sprintf('Time %d',t)); drawnow;

function [ds]=FindDistance(car,target)
dx=target.x-car.x;
dy=target.y-car.y;
ds=sqrt(dx^2+dy^2);

function [dq]=FindAngular(car,target)
dx=target.x-car.x;
dy=target.y-car.y;
q=atan2(dy,dx)-car.q;
dq=AngleWrap(q);

function [q]=AngleWrap(q)
while q<pi,
  q=q+2*pi;
end;
while q>pi,
  q=q-2*pi;
end;

function [car]=CarNew(field)
car.x=(rand/2+0.5)*field.range*sign(randn)*0.8;
car.y=(rand/2+0.5)*field.range*sign(randn)*0.8;
car.q=AngleWrap(randn*2*pi);
car.trace=[car.x; car.y; car.q];
car.shape=[ 2 0; 1 1; -1 1; -1 -1; 1 -1; 2 0]';
car.hdL.shape=plot(car.shape(1,:),car.shape(2,:),'color','b','linewidth',2);
car.hdL.trace=plot(car.trace(1,:),car.trace(2,:),'color',[0 0.66 0]);

function [car]=CarNow(car,time,v,w)
car.x=car.x+time.dt*v*cos(car.q);
car.y=car.y+time.dt*v*sin(car.q);
car.q=car.q+time.dt*w;
car.trace(:,end+1)=[car.x; car.y; car.q];

function [car]=CarShow(fig,car,t)
Rz=[  cos(car.q) -sin(car.q); 
      sin(car.q)  cos(car.q)];
shape=Rz*car.shape+repmat([car.x;car.y],1,6);
set(car.hdL.shape,'xdata',shape(1,:),'ydata',shape(2,:)); 
set(car.hdL.trace,'xdata',car.trace(1,:),'ydata',car.trace(2,:));
title(sprintf('Time %d',t)); drawnow;

function []=CarAxes(car)
miX=min(car.trace(1,:));
mxX=max(car.trace(1,:));
miY=min(car.trace(2,:));
mxY=max(car.trace(2,:));
axis([miX mxX miY mxY]+[-1 1 -1 1]*2); axis equal;
