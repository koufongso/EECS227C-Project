anchorN = 4;
N = 2;
scale = 0.5;
rng(0)
z = (rand(anchorN,2))*scale;
t = sdpvar(N,1);
x = sdpvar(N,2);
h = sdpvar(N,anchorN);
M = 1e9;
cost = M*(t(1)+t(2));

%constraint = [(((z(1,1)-x(1,1))^2+(z(1,2)-x(1,2))^2))*(((z(2,1)-x(1,1))^2+(z(2,2)-x(1,2))^2))*(((z(3,1)-x(1,1))^2+(z(3,2)-x(1,2))^2))<=t(1)];

%  constraint = [(((z(1,1)-x(1,1))^2+(z(1,2)-x(1,2))^2))*(((z(2,1)-x(1,1))^2+(z(2,2)-x(1,2))^2))*(((z(3,1)-x(1,1))^2+(z(3,2)-x(1,2))^2))<=t(1);
%      (((z(1,1)-x(2,1))^2+(z(1,2)-x(2,2))^2))*(((z(2,1)-x(2,1))^2+(z(2,2)-x(2,2))^2))*(((z(3,1)-x(2,1))^2+(z(3,2)-x(2,2))^2))<=t(2)];

% hard code 4 fix, 2 var
 constraint = [(((z(1,1)-x(1,1))^2+(z(1,2)-x(1,2))^2))*(((z(2,1)-x(1,1))^2+(z(2,2)-x(1,2))^2))*(((z(3,1)-x(1,1))^2+(z(3,2)-x(1,2))^2))*(((z(4,1)-x(1,1))^2+(z(4,2)-x(1,2))^2))<=t(1);
     (((z(1,1)-x(2,1))^2+(z(1,2)-x(2,2))^2))*(((z(2,1)-x(2,1))^2+(z(2,2)-x(2,2))^2))*(((z(3,1)-x(2,1))^2+(z(3,2)-x(2,2))^2))*(((z(4,1)-x(2,1))^2+(z(4,2)-x(2,2))^2))<=t(2)];


% for i=1:N
%     offset = (i-1)*N;
%     for j = 1:anchorN
%         h(i,j)=(z(j,1)-x(i,1))^2+(z(j,2)-x(i,2))^2;
%     end
%     constraint = [constraint; prod(h(i,1:anchorN))<=t(i)];
% end
rng(0)
assign(x,z(1:N,:)+(randn(N,2))*0.01,1);
z(1:N,:)
value(x)
options=sdpsettings('solver','ipopt','verbose',0,'showprogress',1,'usex0',1);
result=optimize(constraint,cost,options);
flag=result.problem;
disp(yalmiperror(flag));

xopt = double(x);
%%
figure
hold on
plot(z(:,1),z(:,2),'bo');
plot(xopt(:,1),xopt(:,2),'rx');


%%
t1 = (((z(1,1)-xopt(1,1))^2+(z(1,2)-xopt(1,2))^2))*(((z(2,1)-xopt(1,1))^2+(z(2,2)-xopt(1,2))^2))*(((z(3,1)-xopt(1,1))^2+(z(3,2)-xopt(1,2))^2))*(((z(4,1)-xopt(1,1))^2+(z(4,2)-xopt(1,2))^2))
t2 = (((z(1,1)-xopt(2,1))^2+(z(1,2)-xopt(2,2))^2))*(((z(2,1)-xopt(2,1))^2+(z(2,2)-xopt(2,2))^2))*(((z(3,1)-xopt(2,1))^2+(z(3,2)-xopt(2,2))^2))*(((z(4,1)-xopt(2,1))^2+(z(4,2)-xopt(2,2))^2))

xopt = z(1:N,:);
t1 = (((z(1,1)-xopt(1,1))^2+(z(1,2)-xopt(1,2))^2))*(((z(2,1)-xopt(1,1))^2+(z(2,2)-xopt(1,2))^2))*(((z(3,1)-xopt(1,1))^2+(z(3,2)-xopt(1,2))^2))*(((z(4,1)-xopt(1,1))^2+(z(4,2)-xopt(1,2))^2))
t2 = (((z(1,1)-xopt(2,1))^2+(z(1,2)-xopt(2,2))^2))*(((z(2,1)-xopt(2,1))^2+(z(2,2)-xopt(2,2))^2))*(((z(3,1)-xopt(2,1))^2+(z(3,2)-xopt(2,2))^2))*(((z(4,1)-xopt(2,1))^2+(z(4,2)-xopt(2,2))^2))



%% prepare surface points and their normal vectors
n=2;
s = 2;
p = generateSamplePoints(n,s);
u = zeros(4*n,2);
plot(p(:,1),p(:,2),'-o');
hold on 
for i=1:4*n
    pt = p(i,:);
    u(i,:) = gradInfNorm_2d(pt)';
    plot([pt(1),pt(1)+u(i,1)],[pt(2),pt(2)+u(i,2)],'-')
end


grid on

%%
N = 10;                  % # of stations
P = 4*n;                % # of surface sample points
R = 10;                 % station radius
D = 2*R;                % max distance between stations

% 
xx = visualizeX([1,-3],N,D*ones(N,1));
plot(xx(:,1),xx(:,2),'-x');
%%
x0 = sdpvar(1,2);       % start station
d = sdpvar(N,1);        % distance between adjacent stations
theta = sdpvar(N,1);    % angle between adjacent stations
b = sdpvar(P,N);        % bianary var of checking normal vector
%x = sdpvar(N+1,2);


constraints=[];
for j = 1:P
    for i = 1:N
        constraints = constraints+[b(j,i)==log(exp((x0+d'*[eye(i),zeros(i,N-i);zeros(N-i,N)]*[cos(theta),sin(theta)]-p(j,:))*u(j,:)')+1)];
    end
end

% count = 0;
% M=10;
% for i=1:P
%     for j = 1:N
%         count=count+1;
%         constraints = constraints+[b(i,j) == exp(M*(x(j,:)-p(i,:))*u(i,:)')];
%     end
% end

constraints = constraints + [0<=d<=D];
for i = 1:P
        constraints = constraints+[sum(b(i,:))>=0.7*(N+1)];
end
% 
% for j = 1:N
%         constraints = constraints+[sum(b(:,j))>=1];
% end

cost = ones(1,N)*d;

assign(x0,[1,-3]);
assign(d,2*ones(N,1));
assign(theta,linspace((2*pi)/N,(2*pi),N)');

%assign(x,[0,-2;1,-2;2,-1;2,0;2,1;2,2;1,2;0,2;-1,2;-2,2;-2,0])
options=sdpsettings('solver','ipopt','verbose',0,'showprogress',1,'usex0',1);
result=optimize(constraints,cost,options);
flag=result.problem;
disp(yalmiperror(flag));

%%
x0_opt = double(x0);
bopt = double(b);
dopt = double(d);
thetaopt = double(theta);
x_opt = zeros(N+1,2);
x_opt(1,:) = x0_opt;
for i=1:N
    x_opt(i+1,:) = x0_opt+d'*[eye(i),zeros(i,N-i);zeros(N-i,N)]*[cos(theta),sin(theta)];
end
figure
hold on
plot(p(:,1),p(:,2),'-x');
plot(x_opt(:,1),x_opt(:,2),'-o')
for i=1:N
    text(x_opt(i,1),x_opt(i,2),num2str(i),'fontsize',20);
end

%%
clc;
clear all;
close all;
%% prepare surface points and their normal vectors
m_=5;               % # of points on each side
m = 5*4;            % # of surface sample points
p = generateSamplePoints(m_,10);
p_ = p;
pn = zeros(m,2);
for i=1:m
    pt = p(i,:);
    v = generateUnitNormalVector(pt)';
    pn(i,1:2) = v; 
end



n_ = 5;                            
n = 4*n_;              % # of stations  
R = 2;                 % station radius

x0 = generateSamplePoints(n_,15);

p = reshape(p',[m*2,1]);
pn = reshape(pn',[m*2,1]);
x0 = reshape(x0',[n*2,1]);

figure
hold on
grid on
rectangle('Position',[-5 -5 10 10],'EdgeColor','b')
for i=1:m
    pt = [p((i-1)*2+1),p((i-1)*2+2)];
    pnt = [pn((i-1)*2+1),pn((i-1)*2+2)];
    plot(pt(1),pt(2),'-ob');
    plot([pt(1),pt(1)+pnt(1)],[pt(2),pt(2)+pnt(2)],'--b');
end
temp = reshape(x0,[2,n])';
plot(temp(:,1),temp(:,2),'-rx');hold on

for i=1:n
    xt = [x0((i-1)*2+1),x0((i-1)*2+2)];
    plot(xt(1),xt(2),'-xr');
    viscircles([xt(1),xt(2)],R);
end
%%
tic
maxIteration = 2e4;
stateThreshold = 1e-5;
x = x0+2*(rand()-0.5);
x_pre = inf*ones(2*n,1);
k = 0;
eta = 0.1;
xx = nan(2*n,maxIteration);
while (norm(x-x_pre)>stateThreshold && k<maxIteration)
    [g,g1,g2,g3,g4] = computeGrad(x,p,pn,R,20,k+1);
    x_pre = x;
    x = x-eta/sqrt(k+1)*g/norm(g);
    k = k+1;
    xx(:,k) = x;
end
toc

%%
tic
maxIteration = 1e4;
stateThreshold = 1e-5;

x = x0+0*(rand()-0.5);
x_pre = inf*ones(2*n,1);
k = 0;
eta = 0.1;
xx = nan(2*n,maxIteration);
rng(1);
a = 10;
b = 10;
c = 1;
lambda = [1; a*ones(m,1);b*ones(n,1);c*ones(n,1)];
while (norm(x-x_pre)>stateThreshold && k<maxIteration)
    [g1,g2,g3,g4] = computeGrad_v2(x,p,pn,R,22,k+1);
    g = [g1,g2,g3,g4]*lambda;
    x_pre = x;
    x = x-eta*g/norm(g);
    k = k+1;
    xx(:,k) = x;
end
toc


%%
figure
rectangle('Position',[-5 -5 10 10],'EdgeColor','b')
hold on
xt = xx(:,k);
xt = reshape(xt,[2,n])';
plot(xt(:,1),xt(:,2),'-rx');hold on
viscircles(xt,R*ones(size(xt,1),1));
plot(p_(:,1),p_(:,2),'-bo');
grid on
%%
h = figure;
h.Visible = 'off';
loops = 100;
multiplier = k/loops;
clear M;
M(k) = struct('cdata',[],'colormap',[]);
for i=1:loops
    xt = xx(:,(i-1)*multiplier+1);
    xt = reshape(xt,[2,n])';
    plot(xt(:,1),xt(:,2),'rx');
    hold on;
    plot(p_(:,1),p_(:,2),'-bo');
    drawnow
    M(i) = getframe;
    hold off;
end

h.Visible = 'on';
movie(M,1,20);

%%
vals = [];
for i=1:k
    vals = [vals;objectiveFunciton(xx(:,k),p,pn,R,10)];
end
figure
plot(vals);

%%
[X,Y] = meshgrid(-10:0.1:-5);
surfc(X,Y,-log(max(abs(X),abs(Y))-5));

