% Implicit Scheme 2D Unsteady heat diffusion

tic;
% Defining the domain

domain_size = 80*10^-3;
N = 81;

delta_x = domain_size/(N-1);
delta_y = domain_size/(N-1);
delta_t = 0.1;

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

% Calculate the Fourier's Number
F = alpha*delta_t/delta_x^2;

% Initializing the Error and iteration counter;
time_error = 10;
time_tolerance = 1e-6;
iterations = 0;


% Gauss-Seidel Method
while time_error> time_tolerance

    T_old_time = T;

    residual = 10;
    GS_tolerance = 10^-3;

    while residual > GS_tolerance   
        T_old = T;

        % This loop is for interior nodes
        for i = 2:N-1
            for j = 2:N-1
                T(i,j) = (F*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1))+T_old_time(i,j))/(1 + 4*F);
            end
        end  
        % This loop is for nodes on the left boundary
        for k = 2:N-1
            T(k,1) = (1/(Lambda + h*delta_x))*(h*delta_x*T_ambient + Lambda*T(k,2));
        end

        % Reapplying the boundary conditions (To make the code faster, we can remove the three lines below, but for the sake of my understanding I have kept it
        T(1,:) = T_top;
        T(N,:) = T_Bottom;
        T(:,N) = T_Right;

        % Calculating the Residual
        residual = norm(T-T_old,2);

    end
    iterations = iterations+1; %(Counts the number of iterations in time marching)

    %Calculating the time error
    time_error = norm(T-T_old_time,2);
    % time_error = max(max(abs((T-T_old_time)./T_old_time)))
end

% disp(T) %Final steady state temperature field

elapsed_time = toc;
disp("Total time taken")
disp(elapsed_time)

% Plotting the contour at steady state
figure;
hold on
[X,Y] = meshgrid(0:delta_x:domain_size, 0:delta_y:domain_size); % We get the x,y coordinates
contourf(X, Y, flipud(T), 30);
colormap(jet); colorbar;
title("Temperature Distribution at Steady State using implicit Scheme");
xlabel("X (m)");
ylabel("Y (m)");


