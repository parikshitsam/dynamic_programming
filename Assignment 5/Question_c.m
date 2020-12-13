function invdist = invdist(PI,nz)
    
    Z_old = randi(10,nz,1);
    Z = transpose(Z_old/sum(Z_old));
    
    
    
    z_prime_old = randi(10,nz,1);
    z_prime = transpose(z_prime_old/sum(z_prime_old));
    while Z ~= z_prime
        Z = z_prime;
        z_prime = Z*PI;            
    end 
    invdist = Z;
end
