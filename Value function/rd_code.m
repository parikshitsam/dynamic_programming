%initializing parameters
s=2;
b= 0.99;
a=0.34;
d=0.025;
Z=[1, 0.1];
nk=100;
maxmiter=2000;
pi = [0.9 0.1; 0.1 0.9];
 
%initialize grids
k= zeros(2,nk);
v= zeros(2,nk);
g= zeros(2,nk);
vext= zeros(2,nk,nk);
vnew= zeros(2,nk,nk);
 
%bounds on capital/savings
klb = 0.01;
kub= 5;
 
%defining grig for k
for m = 1:2
    for i = 1:nk
        k(m,i) = klb+(i-1)*(kub-klb)/(nk-1);
    end
end
 %starting value function iteration
 
enditer =0;
iter=0;
 
while (enditer==0)
    iter=iter+1;
    for i=1:nk
        for j=1:nk
            for m=1:2
                if k(j)<= Z(m)*(k(i)^(a)) + (1-d)*k(i)
 
                   vext(m,i,j)= (((Z(m)*(k(i)^(a))+(1-d)*k(i)-k(j))^(1-s))/(1-s)+b*(pi(m,1)*v(m,j)+ pi(m,2)*v(m,j)));
                else
                    vext(m,i,j)= -1/0;
                end
            end
        end
    end
    for m= 1:2
        for i=1:nk
            a = vext(m,i,:);
            b = reshape(a,[1 100]);
            [vnew(m,i),g(m,i)]=max(b);
%              [vnew(m,i),g(m,i)] = max(max(vext(m,i,:)));
        end
    end
    conver=max(max(abs(v-vnew)));
        fprintf( 'Convergence of V (max) = %.8f\n', conver)
 
        if (conver < (10)^(-3) | iter>maxmiter)
           enditer =1;
        end
end