clc
clear
close all

%First, we define the parameters:
beta = 0.99;
alpha = 0.34;
global sigma
sigma = 2;
delta = 0.025;
Z = [1 0.1];
pi = [0.9 0.1; 0.1 0.9];
nz = 2;
iz = invdist(pi,nz);

%Defining the grid size and specifying maximum number of iterations:
nk = 100;
maxiter = 1000;

%initializing grids
%Since we have two states of nature we initialize a 2x100x100 array

k = [zeros(1,nk);zeros(1,nk)];
v = [zeros(1,nk);zeros(1,nk)];
g = [zeros(1,nk);zeros(1,nk)];
vnew = [zeros(1,nk);zeros(1,nk)];

%defining capital bounds
klb = 0.01;
kub = 5;

%defne grid for k. Here we define the grid for both states of nature.
for j= 1:2
    for i = 1:nk
        k(j,i)=klb+(i-1)*(kub-klb)/(nk-1);
    end
end

%define grid for vext. We need to define a multidimensional array for vext

for i = 1:nk
    vext(:,:,i)= [zeros(1,nk);zeros(1,nk)];
end

%start value function iterations
enditer=0;
iter=0;

while(enditer==0)
    iter=iter+1;
    
    for i=1:nk
        for j=1:nk
            for m = 1:2
                if k(j) <= Z(m)*(k(i)^(alpha))+(1-delta)*k(i) 
                    % For the values of k that are feasible
                  c= Z(m)*k(i)^(alpha)+(1-delta)*k(i)-k(j);
%                   vext(m,i,j) = u(c) + beta*pi(m,:)*v(:,j);
                  vext(m,i,j)= u(c)+beta*(iz(1)*v(1,j)+iz(2)*v(2,j));
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
    %plot(k(1:nk),vnew)
    fprintf('Convergence of V(max) = %.8f\n',conver)
    
    if (conver < 10^(-5) | iter>maxiter)
        enditer=1;
    end
end

%finding the policy function
for i = 1:nk
    for m= 1:2
        policy(m,i)=k(g(m,i));
    end
end


disp('')
disp('----Iterations Complete----')

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


function u = u(c)
    global sigma
    u = (c^(1-sigma))/(1-sigma);
end
