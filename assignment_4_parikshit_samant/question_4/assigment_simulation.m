clc
clear

%Assignment 4 Simulation

%Run the pfi

%Variables Given
beta = 0.99;
sigma = 2;
alpha = 0.34;
A = 1;
delta = 0.025;
nk = 200;
maxiter = 10000;
maxiterpol = 10000;


%initialize grids
k = zeros(1,nk);
v = zeros(1,nk);
g = zeros(1,nk);
vext = zeros(nk,nk);
trueg = zeros(1,nk);
truev = zeros(1,nk);
vnew = zeros(1,nk);

% Steady state capital and bounds on grid
kss = (alpha*beta*A/(1-(beta*(1-delta))))^(1/(1-alpha));
klb = 0.1*kss;
kub = 2*kss;

%defne grid for k
for i=1:nk
    k(i) = klb+(i-1)*(kub-klb)/(nk-1);
end
% for i=1:nk
%     trueg(i) = alpha*beta*k(i)^alpha;
%     truev(i) = Atemp + alpha*log(k(i))/(1-alpha*beta);
% end

%Iterate Value function twice


for i = 1:10
    for i=1:nk
        for j=1:nk
            c= max(A*k(i)^(alpha)+(1-delta)*k(i)-k(j),0);
            vext(i,j)= u(c)+beta*v(j);
        end
    end
    
    for i=1:nk
        [v(i),g(i)]=max(vext(i,:));
    end
end

for i=1:nk
    policy(i) = k(g(i));
end


%Policy iteration

enditerpol = 0;
iterpol = 0;
% outer loop: policy

while(enditerpol==0)
    iterpol = iterpol+1;
    %inner loop: 'evaluation' of value function
    iter = 0;
    enditer = 0;
    while( enditer==0)
        iter=iter+1;
        
        for i=1:nk
            c= max(A*k(i)^(alpha)+(1-delta)*k(i)-k(g(i)),0);
            vnew(i) = u(c)+beta*v(g(i));
        end
        conver = max(abs(v-vnew));
        v = vnew;
        fprintf('Convergence of V (max) = %.8f\n',conver)
        if (conver <10^(-5) | iter>maxiter)
            enditer=1;
        end
    end
    
    %Policy Improvement
    for i = 1:nk
        for j = 1:nk
            c = max(A*k(i)^(alpha)+(1-delta)*k(i)-k(j),0);
            vext(i,j) = u(c)+beta*v(j);
        end
    end
    %this gives new policy
    for i=1:nk
        %[v,g] = max(vext(i,:));
        [v(i),g(i)]=max(vext(i,:));
    end
    for i = 1:nk
        policynew(i) = k(g(i));
    end
    
    converpol = max(abs(policy-policynew));
    policy = policynew;
    disp('___________________________________________')
    fprintf('Convergence of policy (max) = %.8f\n',converpol)
    disp('___________________________________________')
    
    if (converpol <10^(-5) |  iterpol > maxiterpol)
            enditerpol=1;
    end
end

% Simulation

%Parameters and Initializing Values and grids
enditersim = 0;
c = zeros(1,nk-1);
k_sim = zeros(1,nk);
k_sim(1) = 0.1*kss;
iter = 1;


while enditersim == 0
    
    k_1_arr = k_sim(iter)* ones(1,nk);
    diff = abs(k_1_arr-k); 
    
    %Find the k closest to the k_sim
    [~, argmin] = min(diff);
    
    %Find k_t+1 for given kt
    k_sim(iter+1)= policy(argmin);
    
    %Find c for given kt and k_t+1
    c(iter) = max(A*k_sim(iter)^(alpha)+(1-delta)*k_sim(iter)-k_sim(iter+1),0);
    
    %Keep a count of loops
    iter = iter +1;
    if (iter==nk)
        enditersim = 1;
    end
    
end
figure(1)
subplot(1,1,1)
plot((1:nk),k_sim,'r',(1:nk-1),c,'b')
legend('k','c')
xlabel('t')
ylabel('c and k')
title('Simulaing kt and ct')


function u = u(c)
    u = (-1/c);
end
      