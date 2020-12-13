clc
clear

%First, we define the parameters:
beta = 0.99;
alpha = 0.34;
global sigma
sigma = 2;
delta = 0.025;
Z = [1 0.1];
pi = [0.9 0.1; 0.1 0.9];

%Defining the grid size and specifying maximum number of iterations:
nk = 100;
maxiter = 2000;

%initializing grids
k = [zeros(1,nk);zeros(1,nk)];
v = [zeros(1,nk);zeros(1,nk)];
g = [zeros(1,nk);zeros(1,nk)];
vnew = [zeros(1,nk);zeros(1,nk)];

%defining capital bounds
klb = 0.01;
kub = 5;

%discretizing for k
for j= 1:2
    for i = 1:nk
        k(j,i)=klb+(i-1)*(kub-klb)/(nk-1);
    end
end

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
                if k(j) <= Z(m)*(k(i))^(alpha)+(1-delta)*k(i)
                   c = Z(m)*k(i)^alpha + (1-delta)*(k(i)) - k(j);
                   vext(m,i,j) = u(c) + beta*pi(m,:)*v(:,j);
                else 
                    vext(m,i,j) = -1/0;
                end
            end
        end
    end
    
    for m = 1:2
        for i=1:nk
            a= vext(m,i,:);
            b = reshape (a, [1 100]);
            [vnew(m,i),g(m,i)] = max(b);
        end
    end
    conver=max(abs(v-vnew));
    v=vnew;
    fprintf('Convergence of V (max) = %.8f\n', conver)
    
    if(conver < 10^(-5) | iter > maxiter)
        enditer=1;
    end
end

for i = 1:nk
    for m= 1:2
        policy(m,i)=k(g(m,i));
    end
end


disp('')
disp('----Iterations Complete----')

mesh(k, Z, v)
function u = u(c)
    global sigma
    u = (c^(1-sigma))/(1-sigma);
end
