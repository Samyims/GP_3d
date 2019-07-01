clear all
close all;

% Load spheretri package
addpath(genpath('spheretri'));

%%

%--------------Off-line Gaussian process regression with known centroid and
%--------------orientation measurements-----------------------%

%---------------------Simulation part-------------------------%

N = 20 ; %Number of points for every face of cube ( number of of measurementts : 6 * N)
z = gen_cube(N) ; %Generation of points from cube
p = (1/2) * ones(3,1) ; %Centroid of cube
R = 0.0001 ; %Noise variance

%-------Generation of noisy measurements
for n = 1 : size(z,2) 
    r_(n) = sqrt((z(1,n) - p(1))^2 + (z(2,n) - p(2))^2 + (z(3,n) - p(3))^2) ;
    theta_(n) = atan2((z(2,n) - p(2)), (z(1,n) - p(1))) ; %Azimuth measurement
    phi_(n) = atan((z(3,n) - p(3))/sqrt(( z(2,n) - p(2))^2 + (z(1,n) - p(1))^2)) ; %Elevation measurement
    Angle_train(n,:) = [theta_(n), phi_(n)] ; 
    c{n} = [ cos(phi_(n)) .* cos(theta_(n)) ; cos(phi_(n)) .* sin(theta_(n)) ; sin(phi_(n))] ;%Orientation vector
    c_(:,n) = c{n} ;
    z_tr(:,n) =  p + r_(n) * c{n} + chol(R)' * randn(3,1) ;
end

%------------------Estimation part---------------------------%

%--Gaussian process parameters
ro = mean(sqrt((z_tr(1,:) - p(1)).^2 + (z_tr(2,:) - p(2)).^2 + (z_tr(3,:) - p(3)).^2)) ;
sigma_r =  0.1 ;
sigma_f  = 2 ;
l = pi/4;
%--Gaussian process kernel
K = @(A, B) kernel_process(A, B, sigma_f, sigma_r, l) ;
%--Angle tests
[basisVertices, ~] = spheretri(20);      % produces evenly spaced points on the sphere
[theta_test, phi_test, ~] = cart2sph(basisVertices(:,1), basisVertices(:,2)...
    , basisVertices(:,3));  % extracts corresponding spherical angles
Nbtest = length(theta_test);
Angle_test = [theta_test, phi_test] ;
Kff = K(Angle_test, Angle_test) ;

% x = ro + chol(Kff)' * randn(Nbtest, 1) ;
% figure, plot(x) 

%---Estimation of output of gaussian process by MC estimate

J_f = [] ;

for m = 1 : length(Angle_train)
   
    Kzf = K(Angle_train(m,:), Angle_test) * K(Angle_test, Angle_test)^-1 ;
    %--Noise regression
    Rf = K(Angle_train(m,:), Angle_train(m,:)) + 0.01 - K(Angle_train(m,:), Angle_test) * ...
        K(Angle_test, Angle_test)^-1 * K(Angle_test, Angle_train(m,:));
    R_{m} = c_(:,m) * Rf * c_(:,m)' + R * eye(3);
    %--Jacobian computation
    J =  c_(:,m) * Kzf ;
    J_f = [ J_f ; J ] ;
    
end

%--Global noise covariance measurements matrix
R_glob = blkdiag(R_{:}) ;
%---Least squares estimation
est =  inv(inv(Kff) + J_f' * inv(R_glob) * J_f) * J_f' *  inv(R_glob)* (z_tr(:) - repmat(p, size(z,2),1));
%--Predicted points
z_pred(:,1) = p(1) + est .* cos(phi_test) .* cos(theta_test) ;
z_pred(:,2) = p(2) + est .* cos(phi_test) .* sin(theta_test) ;
z_pred(:,3) = p(3) + est .* sin(phi_test) ;

%-------------------------Plotting results----------------------%
figure,
plotcube([1 1 1],[ 0  0  0],.8,[1 1 1]), hold all,...
    plot3(z_tr(1,:),z_tr(2,:), z_tr(3,:),'o r ','LineWidth', 3) ;
hold on,
plot3(z_pred(:,1),z_pred(:,2), z_pred(:,3),'x k ','LineWidth', 5) ;
k = boundary(z_pred);
h = trisurf(k, z_pred(:,1), z_pred(:,2), z_pred(:,3),'Facecolor','blue','FaceAlpha',0.1);
view(3)
xlabel('X-axis(m)'), ylabel('Y-axis(m)'), zlabel('Z-axis(m)')
title("Estimation with $l=\ $" + l + ", $\sigma_f =\ $" + sigma_f + "\ and $\sigma_r =\ $" + sigma_r, 'Interpreter','latex')
legend(h,'Estimated shape')
grid minor
%%
% [phi_ord,ind_phi] = sort(phi_test) ;
% [theta_ord,ind_theta] = sort(theta_test) ;
% figure, 
% plot(phi_ord, est(ind_phi),'->'),title('Predicted radius vs elevation test ') ;
% figure,
% plot(theta_ord, est(ind_theta),'->'),title('Predicted radius vs azimuth test') ;
% grid minor

