% Crank Nicolson Scheme 2D Unsteady heat diffusion

% Defining the domain
tic;
domain_size = 80*10^-3;
N = 81;

delta_x = domain_size/(N-1);
delta_y = domain_size/(N-1);
delta_t = 10^-1;

%Defining the initial Temperature of the domain

initial_temperature = 300;
T_top = 423;
T_Bottom = 323;
T_Right = 473;
T_ambient = 573;

T = zeros(N,N);
T(:,:) = initial_temperature;
T(1,:) = T_top;
T(N,:) = T_Bottom;
T(:,N) = T_Right;

% disp(T) % gives us the starting temperature of the domain

% Convection Parameter
h = 250;
Lambda = 5;
alpha = 5*10^-6; % Thermal diffusivity

% Calculate the Modified Fourier's Number (The Fourier's Number here is different from the one used in implicit scheme)
F = alpha*delta_t/(2*delta_x^2);

% Initializing the Error and iteration counter;
error = 10;
tolerance = 10^-3;
iterations = 0;


% Creating the Matrix A and Vector B

while error> tolerance
    T_old  = T;

    for k = 2:N-1
    T(k,1) = (h*delta_x*T_ambient + Lambda*T(k,2))/(Lambda + h*delta_x);
    end

% Initialzing the Matrix A and Vector B

    A = zeros(N*N,N*N);
    B = zeros(N*N,1);

    for i = 1:N
        for j = 1:N
            node = (i-1)*N + j; % gets nodes value from 1 to N*N 
        
            % Points on the Top boundary, Bottom and Right Boundary (Dirichlet Boundary Conditions)
            % I have first calculated the temperature of the points on Left Boundary and then considered it as an Dirichlet condition
            
            if i == 1 || i == N || j == 1 || j == N
                A(node,node) = 1;
                B(node) = T(i,j);   % Because value at boundary is known

            % Interior Points
            else
                A(node, node) = 1+4*F;        % Central point (i, j)
                A(node, node-1) = -F;         % (i, j-1) left neighbor
                A(node, node+1) = -F;         % (i, j+1) right neighbor
    
                A(node, node-N) = -F;         % (i-1, j) top neighbor
                A(node, node+N) = -F;         % (i+1, j) bottom neighbor
                B(node) = (1-4*F)*T(i,j)+ F*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1));  % this comes from the discretized equation
            
            end
        end
    end
    
    T_current = A\B;
    
    Temperature_field = reshape(T_current, [N, N]).';
    % disp(Temperature_field) % gives the temperature field after the first time step

    %Counting the number of iterations
    iterations = iterations + 1;

    %Calculating the temperature difference between two time steps
    Temperature_difference_matrix = Temperature_field-T_old;
    error = norm(Temperature_difference_matrix,2);

    T = Temperature_field;
end
time = toc;
disp(time)
% Plotting the contour at steady state
figure;
hold on
[X,Y] = meshgrid(0:delta_x:domain_size, 0:delta_y:domain_size); % We get the x,y coordinates
contourf(X, Y, flipud(T), 30);
colormap(jet); colorbar;
title("Temperature Distribution at Steady State using Crank-Nicolson Scheme");
xlabel("X (m)");
ylabel("Y (m)");


