clc
clear

%Assignment 4 pfi

%Variables Given
beta = 0.99;
sigma = 2;
alpha = 0.34;
A = 1;
delta = 0.025;
nk = 100;
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
            %Max to get non negative values of C
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

%Plot the Figures
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

% Define Functions
function u = u(c)
    u = (-1/c);
end
        