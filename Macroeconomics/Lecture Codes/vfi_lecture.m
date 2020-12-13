
%Code for VFI

clc
clear
close all
global sigma

%start timer to measure execution time
tic

%parameters

sigma=1;
alpha=0.3;
beta=0.95;
Atemp = 1/(1-beta)*(log(1-alpha*beta)+alpha*beta*log(alpha*beta)/(1-alpha*beta));

%define grid size for k
nk = 100;
maxiter = 5000;

%initialize grids
k = zeros(1,nk);
v = zeros(1,nk);
g = zeros(1,nk);
vext = zeros(nk,nk);
trueg = zeros(1,nk);
truev = zeros(1,nk);
vnew = zeros(1,nk);

% Steady state capital and bounds on grid
kss = (alpha*beta)^(1/(1-alpha));
klb = 0.1*kss;
kub = 2*kss;

%defne grid for k, compute true value and policy
for i=1:nk
    k(i) = klb+(i-1)*(kub-klb)/(nk-1);
end
for i=1:nk
    trueg(i) = alpha*beta*k(i)^alpha;
    truev(i) = Atemp + alpha*log(k(i))/(1-alpha*beta);
end

%Start value function iterations 
enditer = 0;
iter = 0;

while(enditer==0)
    iter=iter+1;
    
    for i=1:nk
        for j=1:nk
            c= max(k(i)^(alpha)-k(j),0);
            vext(i,j)= log(c)+beta*v(j);
        end
    end
    
    for i=1:nk
        [vnew(i),g(i)]=max(vext(i,:));
    end
    conver = max(abs(v-vnew));
    v=vnew;
    fprintf('Convergence of V(max) = %.8f\n',conver)
    
    if (conver < 10^(-5) | iter>maxiter)
        enditer=1;
    end
end

for i=1:nk
    policy(i) = k(g(i));
end

distpol = max(abs(policy-trueg));
disp('')
disp('----- Iterations Complete ---')
fprintf('Distance of policy from true policy (max) = %.8f \n',distpol)
distv = max(abs(v-truev));
fprintf('Distance of V from true V (max) = %.8f\n',distv)
disp(' ')

toc 

figure(1)
subplot(1,2,1)
plot(k(1:nk),trueg(1:nk),'r',k(1:nk),k(1:nk),'k',k(1:nk),policy(1:nk),'b','LineWidth',2)
subplot(1,2,2)
plot(k(1:nk),truev(1:nk),'*r',k(1:nk),v(1:nk),'LineWidth',2)
axis([k(1) k(nk) truev(1) truev(nk)])
           
