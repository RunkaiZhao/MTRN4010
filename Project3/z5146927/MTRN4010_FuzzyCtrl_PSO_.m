function []=MTRN4010_FuzzyCtrl_PSO_()
% fuzzy control
% sim car follow virtual target to target
% PSO optimized initial virtual target location
% objective is to minimise car pose error w.r.t. target

clc; clear all; close all; dbstop if error; warning off;
set(0,'defaultaxesfontname','times new roman');

field.range=50;
time.dt=1; time.T=500; 

[fig]=FigureNew(field);
[car]=CarNew(field); 
[car]=CarNow(car,time,0,0); 
[car]=CarShow(fig,car,0);
[target]=TargetNew(car);
carInit=car; [carInit]=CarShow(fig,carInit,0);
targetInit=target; [targetInit]=TargetShow(fig,target,0);

fis_vel=readfis('MTRN4010_vel.fis');
fis_ang=readfis('MTRN4010_ang.fis');

% define PSO
% you have to decide the number of generations, particles
% the inertia factor, random gain factors
PSO.XLB=5; PSO.XUB=15;
PSO.D=1; PSO.G=10; PSO.N=10;
PSO.w=abs(randn)/2; PSO.cp=2; PSO.cg=2; 
PSO.X=rand(PSO.D,PSO.N); 
PSO.X=PSO.XLB+PSO.X*(PSO.XUB-PSO.XLB);
PSO.Gbest=[]; PSO.gbest=realmax;
PSO.Pbest=PSO.X; PSO.pbest=realmax*ones(PSO.D,PSO.N);
PSO.V=zeros(PSO.D,PSO.N); 

% PSO iterations
for g=1:PSO.G,% generations
  disp(PSO.X);
  for n=1:PSO.N,% particles
    car=carInit; target=targetInit;
    ofs=PSO.X(1,n);% virtual target offset
    % simulate car
    for t=0:time.dt:time.T,
      [target]=TargetNow(t,time,target,ofs);
      [ds]=FindDistance(car,target);
      vel=evalfis(ds,fis_vel); vel=vel*2;
      [dq]=FindAngular(car,target);
      ang=evalfis(dq,fis_ang); ang=ang*2;
      [car]=CarNow(car,time,vel,ang);
    end;
    dx=car.x-target.x; dy=car.y-target.y;
    dq=car.q-target.q; dq=WrapToPi(dq);
    fit(n)=sqrt(dx^2+dy^2+dq^2);
    if fit(n)<PSO.gbest,% update gbest
      PSO.gbest=fit(n); PSO.Gbest=PSO.X(n);
    end;
    if fit(n)<PSO.pbest(n),% update pbest
      PSO.pbest(n)=fit(n); PSO.Pbest=PSO.X(n);
    end;
  end;
  % PSO vel update
  rp=rand(PSO.D,PSO.N); rg=rand(PSO.D,PSO.N);
  PSO.V=PSO.w*PSO.V+rp*PSO.cp.*(PSO.Pbest-PSO.X)+...
    rg*PSO.cg.*(repmat(PSO.Gbest,[1,PSO.N])-PSO.X);
  PSO.X=PSO.X+PSO.V; PSO.X=abs(PSO.X);
  % replace out-of-range particles
  z=find(PSO.X<PSO.XLB);
  PSO.X(z)=PSO.XLB+rand(PSO.D,length(z))*(PSO.XUB-PSO.XLB); 
  z=find(PSO.X>PSO.XUB);
  PSO.X(z)=PSO.XLB+rand(PSO.D,length(z))*(PSO.XUB-PSO.XLB); 
  disp(sprintf('Generation %d Gbest %5.3f gbest %8.6f',g,PSO.Gbest,PSO.gbest));
end;
% show optimized result
car=carInit; target=targetInit;
for t=0:time.dt:time.T,
  [target]=TargetNow(t,time,target,PSO.Gbest);
  [target]=TargetShow(fig,target,t);
  [ds]=FindDistance(car,target);
  vel=evalfis(ds,fis_vel); vel=vel*2;
  [dq]=FindAngular(car,target);
  ang=evalfis(dq,fis_ang); ang=ang*2;
  [car]=CarNow(car,time,vel,ang);
  [car]=CarShow(fig,car,t);
end;
CarAxes(car);


function [fig]=FigureNew(field)
fig.fig=figure('units','normalized','position',[0.3 0.2 0.4 0.5]);
axis(1.5*[-1 1 -1 1]*field.range); hold on; grid on; axis equal;
xlabel('x-direction'); ylabel('y-direction');
fig.ax=axis;

function [target]=TargetNew(car)
target.x=-car.x; target.y=-car.y;
target.q=WrapToPi(rand*2*pi);
target.x_=target.x; target.y_=target.y; 
target.q_=target.q;
target.shape=[ 2 0; 1 1; -1 1; -1 -1; 1 -1; 2 0]';
target.hdL.shape=plot(target.shape(1,:),target.shape(2,:),...
  'color','r','linewidth',2);
target.hdL.shape_=plot(target.shape(1,:),target.shape(2,:),...
  'color','k','linewidth',2);
Rz=[  cos(target.q) -sin(target.q); 
      sin(target.q)  cos(target.q)];
shape=Rz*target.shape+repmat([target.x;target.y],1,6);
set(target.hdL.shape,'xdata',shape(1,:),'ydata',shape(2,:)); 
set(target.hdL.shape_,'xdata',shape(1,:),'ydata',shape(2,:)); 

function [target]=TargetNow(t,time,target,ofs)
ds=50*(1-t/time.T)^3; ds=ds*ofs;
target.x=target.x_-ds*cos(target.q);
target.y=target.y_-ds*sin(target.q);

function [target]=TargetShow(fig,target,t)
Rz=[  cos(target.q) -sin(target.q); 
      sin(target.q)  cos(target.q)];
shape=Rz*target.shape+repmat([target.x;target.y],1,6);
set(target.hdL.shape_,'xdata',shape(1,:),'ydata',shape(2,:)); 
axis(fig.ax); title(sprintf('Time %d',t)); drawnow;

function [ds]=FindDistance(car,target)
dx=target.x-car.x;
dy=target.y-car.y;
ds=sqrt(dx^2+dy^2);

function [dq]=FindAngular(car,target)
dx=target.x-car.x;
dy=target.y-car.y;
q=atan2(dy,dx)-car.q;
dq=WrapToPi(q);

function [q]=WrapToPi(q)
while q<pi,
  q=q+2*pi;
end;
while q>pi,
  q=q-2*pi;
end;

function [car]=CarNew(field)
car.x=(rand-0.5)*field.range*2; 
car.y=(rand-0.5)*field.range*2; 
car.q=WrapToPi(randn*pi);
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
axis(fig.ax); title(sprintf('Time %d',t)); drawnow;

function []=CarAxes(car)
miX=min(car.trace(1,:));
mxX=max(car.trace(1,:));
miY=min(car.trace(2,:));
mxY=max(car.trace(2,:));
axis([miX mxX miY mxY]+[-1 1 -1 1]*2); axis equal;

