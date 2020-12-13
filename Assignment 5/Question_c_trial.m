clc
clear all

P= [0.9 0.1; 0.1 0.9];
% nz = 2;
Z = [1 0];
z_prime = [0.25 0.75];

while Z ~= z_prime
    Z = z_prime;
    z_prime = Z*P;            
end
Z

