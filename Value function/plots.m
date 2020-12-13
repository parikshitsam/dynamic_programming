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
                if k(j) <= Z(m)*(k(i)^(alpha))+(1-delta)*k(i) 
                    % For the values of k that are feasible
                  c= Z(m)*k(i)^(alpha)+(1-delta)*k(i)-k(j);
%                   %vext(m,i,j) = u(c) + beta*pi(m,:)*v(:,j);% This is an
%                   alternate value function which I don't think is right.
                  vext(m,i,j)= u(c)+beta*(id(1)*v(1,j)+id(2)*v(2,j));
                else
                    %If k is not feasible then we set the value of vext to
                    %-inf
                    vext(m,i,j) = -1/0;
                end
                
            end
            
        end
    end

    
    for m= 1:2
        for i=1:nk
            % Since the max function does not work on multidimensional
            % arrays, we reshape the vext to get a 2 dimensional array. 
            
            a = vext(m,i,:);
            b = reshape(a,[1 100]);
            [vnew(m,i),g(m,i)]=max(b);
        end
    end
    
    conver = max(abs(v-vnew));
    v=vnew;
    fprintf('Convergence of V(max) = %.8f\n',conver)
    
    if (conver < 10^(-5) | iter>maxiter)
        enditer=1;
    end
end

%finding the policy function
for i=1:nk
    for m = 1:2
        policy(m,i) = k(g(m,i));
    end
    
end

%Plot for value function
subplot(1,2,1)
mesh(k,Z,v)
xlabel('k')
ylabel('z')
zlabel('v')
title('Value Function')
subplot(1,2,2)

%Plot for policy function
mesh(k,Z,policy)
xlabel('k')
ylabel('z')
zlabel('g(k)')
title('Policy Function')


%Functions
function u = u(c)
    global sigma
    u = (c^(1-sigma))/(1-sigma);
end