function invdist = invdist(PI,nz)
    
    % Initializing the variables
    Z_old = randi(10,nz,1);
    Z = transpose(Z_old/sum(Z_old));
    z_prime_old = randi(10,nz,1);
    z_prime = transpose(z_prime_old/sum(z_prime_old));
    enditer = 0;
    
    while enditer == 0
        Z = z_prime;
        z_prime = Z*PI;
        conver = max(abs(Z - z_prime));
        if conver<0.000005
            enditer = 1;
        end
    end
    
    invdist = Z;
    
end
