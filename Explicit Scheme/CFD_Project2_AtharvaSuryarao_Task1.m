% Explicit Scheme 2D Unsteady heat diffusion

%Defining the domain

domain_size = 80*10^-3;
N = 81;

delta_x = domain_size/(N-1);
delta_y = domain_size/(N-1);
delta_t = 0.05;

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

% Initializing the Error and iteration counter;
error = 10;
tolerance = 10^-3;
iterations = 0;


% Check the stability criteria
F = alpha*delta_t/delta_x^2; % Fouriers Number

if F<=0.25
    disp("Stability satisfied")

    % We use the Gauss Siedell Method as our iterative approach
    T_current = zeros(N,N);
    
    % disp(T_current)
    
    while error>tolerance
        T_old = T;
    
        %This loop is for interior nodes
        for i = 2:N-1
            for j = 2:N-1
                T_current(i,j) = T_old(i,j) + F*(T_old(i+1,j)+T_old(i-1,j)+T_old(i,j+1)+T_old(i,j-1)-4*T_old(i,j));
            end
        end
       
        % This loop is for nodes on Left Boundary
        for k = 2:N-1
            T_current(k,1) = (1/(Lambda + h*delta_x))*(h*delta_x*T_ambient + Lambda*T_old(k,2));
        end
    
        T_current(1,:) = T_top;
        T_current(N,:) = T_Bottom;
        T_current(:,N) = T_Right;
    
    
        %Counting the number of iterations
        iterations = iterations + 1;
    
        %Calculating the Residual
    
        Temperature_difference_matrix = T_current-T_old;
        error = norm(Temperature_difference_matrix,2);
    
        T = T_current;
        % disp(T) % here we can see the evolution of T with time
    end
        
    disp(T)
    
    % Plotting the contour at steady state
    
    figure;
    hold on
    [X,Y] = meshgrid(0:delta_x:domain_size, 0:delta_y:domain_size); % We get the x,y coordinates
    contourf(X, Y, flipud(T), 30);
    colormap(jet); colorbar;
    title("Temperature Distribution at Steady State for Explicit Scheme");
    xlabel("X (m)");
    ylabel("Y (m)");

else
    disp("Stability not satisfied")
end



