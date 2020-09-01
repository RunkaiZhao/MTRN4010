function []=MTRN4010_CarModel_()

clc; clear all; close all; dbstop if error;

% generate/simulate car movement
Rng=30;
T=1000; dt=0.5; 
x=randn*Rng; y=randn*Rng; q=randn*pi;
% prpare graphics
fig1=MyFigure;
[Car.Pose,Car.Shape]=CarShape(x(1),y(1),q(1));
hdL.Pose=plot(x,y);
grid on; hold on; axis equal;
hdL.Shape=plot(Car.Shape.x,Car.Shape.y);
% vel and turn command bounds 
vLB=0; vUB=0.2;
wLB=-pi/20; wUB=+pi/15;
% vel and turn command
[v,w]=GetVW(T,vLB,vUB,wLB,wUB);
% simulation loop
[x,y,q]=MathLoop(dt,T,x,y,q,v,w); %process model
[Car]=ShowPath(hdL,Car,T,Rng,x,y,q);
% collected neural network training data
dx=diff(x); dy=diff(y); dq=WrapToPi(diff(q));
x(end)=[]; y(end)=[]; q(end)=[];
PlotThreeVar(dx,dy,dq,'dx','dy','d\theta','Target data');



% NN model of car movement
inputs=[v; w; q];
targets=[dx; dy; dq];
% network structure, you hvae to decide the number of neurons in the hidden
% layer
hiddenLayerSize = 20; % 
net=fitnet(hiddenLayerSize);
net=configure(net,inputs,targets);
% training arrangement, you have to decide how to arrange training samples
net.divideParam.trainRatio = 0.7;
net.divideParam.valRatio = 0.15;
net.divideParam.testRatio = 0.15;
% train network
[net,tr]=train(net,inputs,targets); % return network and train record
% evaluate network
outputs=net(inputs);
errors=gsubtract(outputs,targets);
performance=perform(net,targets,outputs);
figure; plotperform(tr);
for k=1:size(errors,1),
 figure; ploterrhist(errors(k,:),'bins',30);
end;
disp(net); %view(net);

% supposed real-world loop
nvfig2=MyFigure;
% drive commands, you have to use commands different from traingin
vLB=0; vUB=0.2;
wLB=-pi/20; wUB=+pi/15;
[v,w]=GetVW(T,vLB,vUB,wLB,wUB);
[x,y,q]=MathLoop(dt,T,x,y,q,v,w);
plot(x,y,'b'); hold on; grid on; axis equal;  drawnow;
x_o=x; y_o=y; q_o=q;
for k=1:T-1,% use NN model
  outputs=net([v(k); w(k); q(k)]);
  x(k+1)=outputs(1)+x(k); 
  y(k+1)=outputs(2)+y(k); 
  q(k+1)=outputs(3)+q(k);
  q(k+1)=q(k+1);
  if mod(k,100)==0,
    title(sprintf('%5.0f',k)); drawnow;
  end;
end;
title(sprintf('%5.0f',k+1));
plot(x,y,'r'); drawnow;
hdL=legend('Maths','NN','location','best');
[Car.Pose,Car.Shape]=CarShape(x_o(end),y_o(end),q_o(end));
plot(Car.Shape.x,Car.Shape.y,'b');
[Car.Pose,Car.Shape]=CarShape(x(end),y(end),q(end));
plot(Car.Shape.x,Car.Shape.y,'r'); drawnow;
set(hdL,'string',{'Maths','NN'});
miX=min(x); miY=min(y); mxX=max(x); mxY=max(y);
axis([miX mxX miY mxY]+[-1 1 -1 1]);
xlabel('x-dir'); ylabel('y-dir'); drawnow;
PlotThreeVar(x_o-x,y_o-y,q_o-q,'x-err','y-err','\theta-err','Comparison');


function [x,y,q]=MathLoop(dt,T,x,y,q,v,w)
for k=1:T-1,
  x(k+1)=x(k)+v(k)*dt*cos(q(k));
  y(k+1)=y(k)+v(k)*dt*sin(q(k));
  q(k+1)=q(k)+w(k)*dt;
  q(k+1)=q(k+1);
end;

function [v,w]=GetVW(T,vLB,vUB,wLB,wUB)
v=rand(1,T-1)*(vUB-vLB)+vLB;
w=rand(1,T-1)*(wUB-wLB)+wLB; 

function [Car]=ShowPath(hdL,Car,T,Rng,x,y,q)
for k=1:T-1,
  [Car.Pose(k+1,:),Car.Shape(k+1)]=...
    CarShape(x(k+1),y(k+1),q(k+1));
  set(hdL.Pose,'xdata',Car.Pose(:,1),'ydata',Car.Pose(:,2));
  set(hdL.Shape,'xdata',Car.Shape(k+1).x,'ydata',Car.Shape(k+1).y);
  if mod(k,100)==0,
    title(sprintf('%5.0f',k));
  end;
  drawnow;
end;
title(sprintf('%5.0f',k+1));
miX=min(Car.Pose(:,1)); miY=min(Car.Pose(:,2));
mxX=max(Car.Pose(:,1)); mxY=max(Car.Pose(:,2));
plot(Car.Pose(:,1),Car.Pose(:,2),'b')
axis([miX mxX miY mxY]+[-1 1 -1 1]); axis equal;
xlabel('x-dir'); ylabel('y-dir'); drawnow;

function []=PlotThreeVar(dx,dy,dq,txt_x,txt_y,txt_q,txt_title);
figure; 
subplot(3,1,1); plot(dx); grid on; ylabel(txt_x); title(txt_title);
subplot(3,1,2); plot(dx); grid on; ylabel(txt_y);
subplot(3,1,3); plot(dq); grid on; ylabel(txt_q); xlabel('Sample');
drawnow;

function [miX,mxX,miY,mxY]=FinalPlot(Car)
miX=min(Car.Pose(:,1)); miY=min(Car.Pose(:,2));
mxX=max(Car.Pose(:,1)); mxY=max(Car.Pose(:,2));
plot(Car.Pose(:,1),Car.Pose(:,2),'b')
axis([miX mxX miY mxY]+[-1 1 -1 1]); axis equal;
xlabel('x-dir'); ylabel('y-dir');

function [CarPose,CarShape]=CarShape(x,y,q)
L=1; W=0.5;
CarPose=[x y q];
Rot=[cos(q) -sin(q); sin(q) cos(q)];
xy=Rot*[L 0 0 L; 0 W -W 0];
CarShape.x=xy(1,:)+x;
CarShape.y=xy(2,:)+y;

function [fig]=MyFigure()
fig=figure;
set(fig,'units','normalized','position',...
  [rand*0.2+0.2 rand*0.2+0.2 0.30 0.45]);


function [q]=WrapToPi(q)
while q<pi,
  q=q+2*pi;
end;
while q>pi,
  q=q-2*pi;
end;
