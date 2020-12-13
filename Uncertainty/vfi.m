clc
clear
close all

global sigma
sigma = 2;


%Variables

beta = 0.99;
alpha = 0.34;
Z = [1 0.1];
delta = 0.025;
nk = 100;
maxiter = 10000; %The function converges in 1233 iterations


%initialize grids
%Since we have two states of nature we initialize a 2x100 array

k = zeros(2,nk);
v = zeros(2,nk);
g = zeros(2,nk);
vext = zeros(2,nk,nk);
vnew = zeros(2,nk);
pi = [0.9 0.1; 0.1 0.9];
nz = 2;
id = invdist(pi,nz);

% Upper and Lower Bounds
klb = 0.01;
kub = 5;

%defne grid for k. 
%%Here we define the grid for both states of nature.

for m = 1:2
    for i=1:nk
        k(m, i) = klb+(i-1)*(kub-klb)/(nk-1);
    end
end


%Start value function iterations 
enditer = 0;
iter = 0;

while(enditer==0)
    iter=iter+1;
    
    for i=1:nk
        for j=1:nk
            for m = 1:2
            c= max(A*k(i)^(alpha)+(1-delta)*k(i)-k(j),0);
            vext(i,j)= u(c)+beta*v(j);
        end
    end
    
    for i=1:nk
        [vnew(i),g(i)]=max(vext(i,:));
    end
    conver = max(abs(v-vnew));
    v=vnew;
    %plot(k(1:nk),vnew)
    fprintf('Convergence of V(max) = %.8f\n',conver)
    
    if (conver < 10^(-5) | iter>maxiter)
        enditer=1;
    end
end

for i=1:nk
    policy(i) = k(g(i));
end

%Plot figures
figure(1)
subplot(1,2,1)
plot(k(1:nk),k(1:nk),'k',k(1:nk),policy(1:nk),'b','LineWidth',2)
legend('k',"approximate k'")
title('Policy Function')
xlabel('k') 
ylabel('g(k)')
subplot(1,2,2)
plot(k(1:nk),v(1:nk),'LineWidth',2)
legend('approximate v')
title('Value Function')
xlabel('k') 
ylabel('v(k)')


%Functions
function u = u(c)
    u = (c^(1-2))/(1-2);
end
