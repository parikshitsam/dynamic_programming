clc
clear
tic

beta = 0.99;
alpha = 0.34;
global sigma
sigma = 2;
delta = 0.025;
Z = [1 0.1];
pi = [0.9 0.1; 0.1 0.9];



nk = 100;
maxiter = 5000;

%initializing grids
k = zeros(1,nk);
v = zeros(2,nk);
g = zeros(2,nk);
vnew = zeros(2,nk);


klb = 0.01;
kub = 5;

for i=1:nk
    k(i)=klb+(i-1)*(kub-klb)/(nk-1);
end

% Start Value Function Iterations
enditer=0;
iter=0;

while(enditer==0)
    iter=iter+1;
    
    for m = 1:2
        for i=1:nk
            for j=1:nk
                c= max(Z(m)*(k(i))^alpha + (1-delta)*(k(i)) - k(j),0);
                vext(m,i,j) = u(c) + beta*pi(m,:)*v(:,j);
            end
        end
    end
    
    for m=1:2
        for i=1:nk
        [vnew(m,i),g(m,i)]=max(vext(m,i,:));
        end
    end
    
    conver=max(abs(v-vnew));
    v=vnew;
    fprintf('Convergence of V (max) = %.8f\n', conver)
    
    if(conver < 10^(-5) | iter > maxiter)
        enditer=1;
    end
end


for m=1:2
    for i=1:nk
        policy(m,i)=k(g(m,i));
    end
end

disp('')
disp('----Iterations Complete----')

toc

subplot(1,2,1)
mesh(k,Z,v)

subplot(1,2,2)
mesh(k,Z,policy)
%Functions
function u = u(c)
    global sigma
    u = (c^(1-sigma))/(1-sigma);
end