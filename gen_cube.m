function z = gen_cube(N)
%--Function generates points from cube drawn
% every face with an uniform distribution
x_1  = rand(1, N) ;
y_1  = rand(1, N) ;
z_1  = zeros(1, N) ;

x_2  = zeros(1, N) ;
y_2  = rand(1, N) ;
z_2  = rand(1, N) ;

x_3 = rand(1, N) ;
y_3 = zeros(1, N) ;
z_3 = rand(1, N) ;

x_4  = rand(1, N) ;
y_4  = rand(1, N) ;
z_4  = ones(1, N) ;

x_5 = rand(1, N) ;
y_5 = ones(1, N) ;
z_5 = rand(1, N) ;

x_6  = ones(1, N) ;
y_6  = rand(1, N) ;
z_6  = rand(1, N) ;

z = cat( 2,[[ x_1 ; y_1 ; z_1 ], [x_2 ; y_2 ; z_2], ...
    [x_3 ; y_3 ; z_3] , [x_4 ; y_4 ; z_4], ...
    [x_5 ; y_5 ; z_5] , [x_6 ; y_6 ; z_6] ] );

end

