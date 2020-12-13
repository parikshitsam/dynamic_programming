clc
clear
close all

global sigma
sigma = 2;


%Variables
%Variables Given
beta = 0.99;

alpha = 0.34;
Z = [1 0.1];
delta = 0.025;
nk = 100;
maxiter = 10000;

%initialize grids
k = [zeros(1,nk);zeros(1,nk)];
v = [zeros(1,nk);zeros(1,nk)];
g = [zeros(1,nk);zeros(1,nk)];
vext = [zeros(1,nk);zeros(1,nk)];
vnew = [zeros(1,nk);zeros(1,nk)];
iz = [0.5 0.5];

% Steady state capital and bounds on grid
klb = 0.01;
kub = 5;

%defne grid for k
for m = 1:2
    for i=1:nk
        k(m, i) = klb+(i-1)*(kub-klb)/(nk-1);
    end
end
%define grid for vext
for i = 1:nk 
    vext(:,:,i) = [zeros(1,nk);zeros(1,nk)];
end

% vext(2,1,1)
%Start value function iterations 
enditer = 0;
iter = 0;

while(enditer==0)
    iter=iter+1;
    
    for i=1:nk
        for j=1:nk
            for m = 1:2
                if k(j) <= Z(m)*k(i)^(alpha)+(1-delta)*k(i)
                  c= Z(m)*k(i)^(alpha)+(1-delta)*k(i)-k(j);
                  vext(m,i,j)= u(c)+beta*(iz(1)*v(m,j)+iz(2)*v(m,j));
                else
                    vext(m,i,j) = -1/0;
                end
                
            end
            
        end
    end
    
%     a=vext(1,1,:);
%     b = reshape(a,[1 100]);
    for m= 1:2
        for i=1:nk
            a = vext(m,i,:);
            b = reshape(a,[1 100]);
            [vnew(m,i),g(m,i)]=max(b);
        end
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
    for m = 1:2
        policy(m,i) = k(g(m,i));
    end
    
end

mesh(k,Z,v)
% mesh(policy,Z,k)

%Functions
function u = u(c)
    global sigma
    u = (c^(1-sigma))/(1-sigma);
end