function K = kernel_process( Gam_train, Gam_test, sigma_f, sigma_r, l )
%--Function computes the GP kernel for tow set of inputs
%Gam_test and Gam_train with parameters sigma_f, sigma_r and l
for i = 1 : size(Gam_train, 1)
    
    for j  = 1 : size(Gam_test, 1)
        
        K(i,j) = kernel(Gam_train(i,:),Gam_test(j,:), sigma_f, sigma_r, l) ;
        
    end
    
end


end

function k = kernel(Gam1, Gam2, sigma_f, sigma_r, l)

c = @(Gam1, Gam2) (acos(cos(Gam1(2))*cos(Gam2(2))*cos(Gam1(1))*cos(Gam2(1))+ ...
cos(Gam2(2))*cos(Gam1(2))*sin(Gam1(1))*sin(Gam2(1))+sin(Gam1(2))*sin(Gam2(2))))^2 ;
k = sigma_f^2 * exp(-c(Gam1, Gam2)/(2*l^2)) + sigma_r^2 ;


end