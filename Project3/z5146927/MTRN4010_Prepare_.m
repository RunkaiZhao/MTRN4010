function []=MTRN4010_Prepare_()
% simulate car in random movement

clc; clear all; close all; dbstop if error;
set(0,'defaultaxesfontname','times new roman');

field.range=50;
time.dt=1; time.T=500;

[fig]=FigureNew(field);
[car]=CarNew();

for t=0:time.dt:time.T,
  [car]=CarNow(car,time,rand,randn*0.25);
  [car]=CarShow(car,t);
end;
CarAxes(car);


function [fig]=FigureNew(field)
fig=figure('units','normalized','position',[0.1 0.2 0.4 0.5]);
% axis([-1 1 -1 1]*field.range*1.5); 
hold on; grid on; axis equal;
xlabel('x-direction'); ylabel('y-direction');

function [car]=CarNew()
car.x=0; car.y=0; car.q=0;
car.trace=[car.x; car.y; car.q];
car.shape=[ 2 0; 1 1; -1 1; -1 -1; 1 -1; 2 0]';
car.hdL.shape=plot(car.shape(1,:),car.shape(2,:),'color','b','linewidth',2);
car.hdL.trace=plot(car.trace(1,:),car.trace(2,:),'color',[0 0.66 0]);

function [car]=CarNow(car,time,v,w)
% you have to imlement the mathematical model fro the car
% and wrap the car orientation usig the WrapToPi function
car.trace(:,end+1)=[car.x; car.y; car.q];

function [car]=CarShow(car,t)
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

function [q]=WrapToPi(q)
while q<pi,
  q=q+2*pi;
end;
while q>pi,
  q=q-2*pi;
end;


