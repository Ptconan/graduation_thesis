%%%%% This code is to optimize the truss structure design using SDP,
%%%%% complemented by strmM.txt and sqrtK.txt
%%%%% By Mengmeng Song, 2026/05/19
%%%%% objective: minimizing volumn of the structure/ weight 
%%%%% variable: x, cross sectional areas of bars
%%%%% constriants: equilibrium equation, nonegativeness of bars, and 
nd=46;%number of freedom
n=200;%the dimension of variable x, i.e., number of bars in the ground structure
d_nd=2*nd+1;
load('sqrtK.txt');
load('strmM.txt');
gamma = 1;
lambda = 5e6;
fl = zeros(nd,1);
fl(41) = 4;
full_nd=2*nd+1;
tic;
cvx_begin sdp
variable x(n, 1)
dual variables Y mu;
expression Fullmatrix(d_nd)
minimize(trace(ones(1,n) * x))
subject to
Y:[gamma -fl' zeros(1,full_nd-nd-1); -fl sqrtK*sparse(diag(x))*sqrtK' zeros(nd,full_nd-nd-1); zeros(full_nd-nd-1,nd+1) sqrtK*sparse(diag(x))*sqrtK'-lambda*reshape(strmM*x,nd,nd)]>=0;
mu:x>=0;
cvx_end
toc;

[dll,matH,coord_x,ir,irr,ird,Idx] = member(4,4,Inf);
draw_cs_freq(coord_x,irr,x,1,4,4);


A0=[1 -fl' zeros(1,nd); -fl zeros(nd,nd) zeros(nd,nd); zeros(nd,1) zeros(nd,nd) zeros(nd,nd)];
save('myData.mat', 'A0');
for i=1:200
A=[zeros(1,2*nd+1); zeros(nd,1) sqrtK(:,i)*sqrtK(:,i)' zeros(nd,nd); zeros(nd,nd+1) sqrtK(:,i)*sqrtK(:,i)'-lambda*reshape(strmM(:,i),nd,nd)];
expr=['A' num2str(i) '=' 'A'];
eval(expr);
expr=['A' num2str(i)];
save('myData.mat', expr,'-append');
end

load('myData.mat');
AI=zeros(2*nd+1,n*(2*nd+1));
for i=1:n
expr=['A' num2str(i)];
AI(:,(2*nd+1)*(i-1)+1:(2*nd+1)*i)=eval(expr);
end
tic;
cvx_begin sdp
variable x(n, 1)
expression Fullmatrix(full_nd)
Fullmatrix = A0+AI*kron(x,eye(2*nd+1));
Fullmatrix = sparse(Fullmatrix);
minimize(trace(ones(1,n) * x))
subject to
Fullmatrix >= 0;
x >= 0;
cvx_end
toc;