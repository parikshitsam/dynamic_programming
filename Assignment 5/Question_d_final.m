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
pi = [0.9 0.1; 0.1 0.9];
nz = 2;
% iz = invdist(pi,nz);

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
               
                  c= max(Z(m)*k(i)^(alpha)+(1-delta)*k(i)-k(j),0.000001);
                  bt = 0;
                  for p = 1:2
                      bt = bt + (pi(m,p)*v(m,j)+pi(m,p)*v(m,j));
%                       vext(m,i,j)= u(c)+beta*(pi(i,j)*v(m,j)+Z(2)*v(m,j));
                  end
                  vext(m, i , j) = (c^(1-2))/(1-2)+ beta*bt;

                
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
%              [vnew(m,i),g(m,i)] = max(max(vext(m,i,:)));
        end
    end
    
    conver = max(abs(v-vnew));
    v = vnew;
    fprintf('Convergence of V (max) = %.8f\n',conver)
    if (conver <10^(-5) | iter>maxiter)
            enditer=1;
    end
%     conver_n = [0 0];
%     for m = 1:2
%     conver_n(m) = max(abs(v(m)-vnew(m)));
%     end
%     conver = max(conver_n);
%     v=vnew;
%     %plot(k(1:nk),vnew)
%     fprintf('Convergence of V(max) = %.8f\n',conver)
%     
%     if (conver < 10^(-5) | iter>maxiter)
%         enditer=1;
%     end
end

for i=1:nk
    for m = 1:2
        policy(m,i) = k(g(m,i));
    end    
end

mesh(k,Z,v)
% xlabel('k')
% ylabel('z')
% % zlabel('v')
% mesh(policy,Z,k)

%Functions
% function u = u(c)
%     global sigma
%     u = (c^(1-sigma))/(1-sigma);
% end