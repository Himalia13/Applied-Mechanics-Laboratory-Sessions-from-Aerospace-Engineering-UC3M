%Sofıa Marın Puche
%Alejandra Menendez Lucini
%Andres Velazquez Vela


% set the size of the ticks on the axes
set(groot, 'defaultLegendFontSize', 12);
% set the default size of the text
set(groot, 'defaultTextFontSize', 12);
% set the default axes font size.
set(groot, 'defaultAxesFontSize', 12);

% set the width of the axes
set(groot, 'defaultAxesLineWidth', 1);
% activate the minor ticks of the axes
set(groot, 'defaultAxesXMinorTick', 'on');
set(groot, 'defaultAxesYMinorTick', 'on');
% deactivate the legend by default
set(groot, 'defaultLegendBox', 'off');
% define the default line width in the plots
set(groot, 'defaultLineLineWidth', 1);
% define the default line marker size
set(groot, 'defaultLineMarkerSize', 5);
% set the font of the axes ticks to Latex
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
% define the font for the default text for the rest of
% objects (labels, titles, legend, etc...)
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

close all
clear all
clc


g = 9.81; % m/s`2
m = 1; % kg
a = 0.2; % m
xi0 = 1; % m
omega = [0.1, 0, 0, 0, 0.1];
velocityp = [0, -1, -10, 4, 1];
%index = input('Enter the case to analise:');

for index = 1:5
    phi0 = 0; %initial angle
    dphi0 = calcDphi0(a, xi0, omega, velocityp, index); %initial angle velocity
   
    
    %integration span and the initial conditions
    tspan = [0 60];
    X0 = [phi0; dphi0];
    
    %solving the ODE with ode45
    Opt = odeset('Events', @(t,X)stopfun(t,X, xi0, m, g, phi0, omega, a, index), 'Stats', 'on', 'Refine', 20, 'RelTol', 1e-8,'AbsTol', 1e-10);%'Refine', 10
    [t, X] = ode45(@(t, X) diffeq(t, X, g, omega, a, xi0, phi0, index), tspan, X0, Opt);
   
   
    %Plot the evolution of φ and the amount of unspooled string ξ as a function of time
    xi = xi0-omega(index)*t*a + (X(:,1)- phi0)*a;
    
    figure(1+6*(index-1))
    hold on
    
    title(sprintf('Case %d: $\\xi$ in function of time', index), 'Interpreter', 'latex');
    xlabel('time (t)');
    
    h1 = plot(t, xi, 'k');
    set(h1, 'Color', [0.85, 0.33, 0.10]);  
    set(gca, 'YColor', [0.85, 0.33, 0.10]); 
    ylabel('$\xi(m)$');
    legend('$\xi(m)$');
    
    switch index
        case 1
            xlim([0, max(t)]);
        case 2
            ylim([0.8, 1.2]);
            xlim([0, 20]);
        case 3
            ylim([0, 1.1])
            xlim([0, max(t)]);
        case 4
            ylim([0.5, 1.5]);
            xlim([0, max(t)]);
        case 5
            ylim([0, 1.2]);
            xlim([0, max(t)]);
            
    end
    
    grid on;
    %print(1+6*(index-1), '-dpng', '-r600');
    figure(2+6*(index-1))
    hold on
    
    title(sprintf('Case %d: $\\phi$ in function of time.', index), 'Interpreter', 'latex');
    xlabel('time (t)');
    
    
    h2 = plot(t, X(:,1), 'k');
    set(h2,'Color', [0, 0.45, 0.74]);  
    set(gca, 'YColor', [0, 0.45, 0.74]);  
    ylabel('$\phi(rad)$');
    legend('$\phi(rad)$');
    switch index
        case 1
            xlim([0, max(t)]);
        case 2
            ylim([-1, 1]);
            xlim([0, 20]);
        case 3
            ylim([-5, 0]);
            xlim([0, max(t)]);
        case 4
            ylim([-2, 2]);
            xlim([0, max(t)]);
        case 5
            ylim([-1.5, 1]);
            xlim([0, max(t)+0.3]);
    end
    
    grid on;
    %print(2+6*(index-1), '-dpng', '-r600');
    %Plot the tension in the string T as a function of time.
    T = m*(g*cos(X(:,1)) + X(:,2).^2 .* (xi0 - a*(omega(index)*t - (X(:,1) - phi0))));
    
    figure(3+6*(index-1))
    hold on
    h3 = plot(t, T, 'k');
    set(h3, 'Color', [34 / 255, 139 / 255, 34 / 255]); 
    
    switch index
        case 2
            ylim([8, 12]);
            xlim([0, 20]);
        case 3
            ylim([0, 3000]);
            xlim([0, max(t)]);
        case 4
            ylim([0, 30]);
            xlim([0, max(t)]);
        case 5
            ylim([2, 34]);
            xlim([0, max(t)+0.3]);
        case 1
            xlim([0, max(t)]);
    end
    
    grid on;
    
    title(sprintf('Case %d: $T$ in function of time.', index), 'Interpreter', 'latex');
  
    xlabel('time (t)');
    
    legend('$T(N)$');
    
    set(gca, 'YColor', [34 / 255, 139 / 255, 34 / 255]);  
    ylabel('$T(N)$');
    %print(3+6*(index-1), '-dpng', '-r600');

    %Plot the Energy as a function of time.
    V_squared = (X(:,2).*(xi0- a*(omega(index)*t - (X(:,1)- phi0)))).^2 + (a*omega(index)).^2;
    Energy = 0.5*m*V_squared + m*g*(1+ a * sin(X(:,1)) + (-xi0 + a * (omega(index) * t - (X(:,1) - phi0))) .* (cos(X(:,1))));
    
    figure(4+6*(index-1))
    hold on
    h4 = plot(t, Energy, 'k');
    set(h4, 'Color', [139 / 255, 0 / 255, 0 / 255]);

    hold on; 
     h5 = plot(t, m*g*(1+ a * sin(X(:,1)) + (-xi0 + a * (omega(index) * t - (X(:,1) - phi0))) .* (cos(X(:,1)))), 'k--');
    set(h5, 'Color', [34 / 255, 139 / 255, 34 / 255]);
    
     hold on; 
     h6 = plot(t, 0.5*m*V_squared, 'k--');
    set(h6, 'Color', [1, 0, 0]);

    legend('$E(J)$', 'Potential', 'Kinetic');
    switch index
        case 2
            ylim([-0.1, 0.7]);
            xlim([0, 20]);
        case 3
            ylim([0, 51]);
            xlim([0, max(t)]);

            legend({'$E(J)$', 'Potential', 'Kinetic'}, 'Location', 'east');
        case 4
            ylim([-1, 11]);
            xlim([0, max(t)]);
        case 5
            ylim([0, 8]);
            xlim([0, max(t)]);

            legend({'$E(J)$', 'Potential', 'Kinetic'}, 'Location', 'west');
        case 1
            xlim([0, max(t)]);
            ylim([-1, 10]);

            legend({'$E(J)$', 'Potential', 'Kinetic'}, 'Location', 'east');
    end
    
    grid on;
    
    title(sprintf('Case %d: $Energy$ in function of time respect to the point P(-1,0)', index), 'Interpreter', 'latex');
  
    xlabel('time (t)');
    
    
    
    set(gca, 'YColor', [139 / 255, 0 / 255, 0 / 255]);  
    ylabel('$E(J)$');
    %print(4+6*(index-1), '-dpng', '-r600');
    
    % Plotting the trajectory with color gradient
    x = a * cos(X(:,1)) + (-xi0 + a * (omega(index) * t - (X(:,1) - phi0))) .* (-sin(X(:,1)));
    y = a * sin(X(:,1)) + (-xi0 + a * (omega(index) * t - (X(:,1) - phi0))) .* (cos(X(:,1)));

    theta = linspace(0, 2 * pi, 200); % 200 points around the circle
    circle_x = a*cos(theta);           % x-coordinates of the circle
    circle_y = a*sin(theta);           % y-coordinates of the circle

    % Fine time vector for interpolation
    t_fine = linspace(min(t), max(t), 5000); % More points for smooth curve

    % Interpolating x and y using interp1
    [t_unique, ia] = unique(t);
    x_unique = x(ia);
    y_unique = y(ia);
    x_interp = interp1(t_unique, x_unique, t_fine, 'spline'); % 'spline' for smoothness
    y_interp = interp1(t_unique, y_unique, t_fine, 'spline');

    % color gradient based on time
    colorg = linspace(0, max(t), length(t_fine)); % Gradient for interpolated points

    figure(5+6*(index-1));
    hold on;
    scatter(x_interp, y_interp, 16, colorg, 'filled');  
    colormap(jet);  

    xlabel('$x(t)$ (meters)');
    ylabel('$y(t)$ (meters)');
    title(sprintf('Case %d: 2D Trajectory with Time Gradient', index), 'Interpreter', 'latex');


    plot(circle_x, circle_y, 'k', 'LineWidth', 1.5); 
    legend({'Trajectory', 'Unit Circle'}, 'Interpreter', 'latex'); 
    grid on;
    axis equal;


    cb = colorbar;

    cb.Label.FontSize = 14; 
    cb.Ticks = linspace(0, max(t), 5);  
    cb.TickLabels = arrayfun(@(v) sprintf('%.1f', v), cb.Ticks, 'UniformOutput', false);

    cb.Label.String = 'Time (s) in a colorbar'; 
    cb.Label.Interpreter = 'latex'; 
    %print(5+6*(index-1), '-dpng', '-r600'); 
    

    % % Uncomment this section if you want to save th videos. 
    % % Initialize video writer object
    % video_name = sprintf('trajectory_video_case_%d.mp4', index); 
    % v = VideoWriter(video_name, 'MPEG-4');
    % switch index
    %     case 1
    %         v.FrameRate = 100;  % frame rate
    % 
    %         open(v);  % Open video file for writing
    % 
    %         % Fine time vector for interpolation
    %         t_fine = linspace(min(t), max(t), 20*max(t)); % More points for smooth curve
    %     case 2
    %         v.FrameRate = 100;  %frame rate
    % 
    %         open(v);  % Open video file for writing
    % 
    %         % Fine time vector for interpolation
    %         t_fine = linspace(min(t), max(t), 100*max(t)); 
    %     case 3
    %         v.FrameRate = 100;  % rame rate
    % 
    %         open(v);  % Open video file for writing
    % 
    %         % Fine time vector for interpolation
    %         t_fine = linspace(min(t), max(t), 1000*max(t)); 
    %     case 4
    %         v.FrameRate = 100;  % frame rate
    % 
    %         open(v);  % Open video file for writing
    % 
    %         % Fine time vector for interpolation
    %         t_fine = linspace(min(t), max(t), 200*max(t)); 
    %     case 5
    %         v.FrameRate = 60;  % frame rate
    % 
    %         open(v);  % Open video file for writing
    % 
    %         % Fine time vector for interpolation
    %         t_fine = linspace(min(t), max(t), 60*max(t)); 
    % end
    % 
    % 
    % 
    % 
    % 
    % % Interpolating x and y using interp1
    % [t_unique, ia] = unique(t);
    % x_unique = x(ia);
    % y_unique = y(ia);
    % x_interp = interp1(t_unique, x_unique, t_fine, 'spline');
    % y_interp = interp1(t_unique, y_unique, t_fine, 'spline');
    % [t_unique, ia] = unique(t);
    % X_unique = X(ia, 1);
    % phi_interp = interp1(t_unique, X_unique, t_fine, 'spline');
    % 
    % num_colors = length(t_fine); 
    % 
    % dblue = [0, 0, 1];
    % blue = [0, 0, 1];    
    % cyan = [0, 1, 1];    
    % green = [0, 1, 0];   
    % yellow = [1, 1, 0];  
    % red = [1, 0, 0];     
    % dred = [0.5, 0, 0];
    % 
    % 
    % jet_colors = [dblue; blue; cyan; green; yellow; red; dred];
    % 
    % 
    % cmap_jet = interp1(linspace(0, 1, size(jet_colors, 1)), jet_colors, linspace(0, 1, num_colors));
    % 
    % total_time = max(t);
    % 
    % 
    % figure_handle = figure(6 + 6 * (index - 1)); 
    % set(figure_handle, 'Position', [0, 0, 1920, 1080]);  
    % grid on;
    % 
    % for i = 1:length(t_fine)
    %     clf;  % Clear figure for each frame
    % 
    %     % mask for the colors based on the current index
    %     alpha_values = zeros(length(t_fine), 1);
    %     alpha_values(1:i) = 1; 
    %     alpha_values(i+1:end) = 0; 
    % 
    % 
    %     scatter(x_interp(1:i), y_interp(1:i), 3, cmap_jet(1:i, :), 'filled');
    %     hold on
    %      scatter(x_interp(i), y_interp(i), 34, cmap_jet(i, :), 'filled');
    % 
    % 
    %     % Plot the unit circle
    %     hold on;
    %     theta = linspace(0, 2*pi, 200);
    %     circle_x = a * cos(theta);
    %     circle_y = a * sin(theta);
    %     plot(circle_x, circle_y, 'k', 'LineWidth', 1.5);  
    % 
    % 
    %     x_tangent = a * cos(phi_interp(1:i)); 
    %     y_tangent = a * sin(phi_interp(1:i));
    % 
    % 
    % 
    %     % Plot the tangent line
    %     plot([x_tangent(1:i), x_interp(i)], [y_tangent(1:i), y_interp(i)], 'k', 'LineWidth', 1);
    % 
    % 
    %     % Axis labels and title
    %     xlabel('$x(t)$ (meters)', 'Interpreter', 'latex');
    %     ylabel('$y(t)$ (meters)', 'Interpreter', 'latex');
    %     title(sprintf('Case %d: 2D Video', index), 'Interpreter', 'latex');
    % 
    % 
    %     axis equal;
    % 
    % 
    %      elapsed_time = t_fine(i); % Get the elapsed time
    %      text(0, 0, sprintf('%.1f s', elapsed_time), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');
    % 
    %      switch index
    %         case 1
    %             text(-0.3, -0.4, sprintf('speed x 5'), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex'); 
    %         case 2
    %            text(-0.3, -0.4, sprintf('speed x 1'), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex'); 
    % 
    %         case 3
    %             text(-0.3, -0.4, sprintf('speed x 0.1'), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex'); 
    %         case 4
    %            text(-0.3, -0.4, sprintf('speed x 0.5'), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');
    %         case 5
    %          text(-0.3, -0.4, sprintf('speed x 1 '), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex'); 
    %      end
    % 
    %     % Capture the frame
    %     frame = getframe(gcf);
    % 
    %     % Write the frame to the video
    %     writeVideo(v, frame);
    % 
    %     if i > 1
    %         pause(t_fine(i) - t_fine(i - 1)); 
    %     end
    % end
    % for i = 1:500
    % 
    %     frame = getframe(gcf);
    % 
    % 
    %     writeVideo(v, frame);
    % end
    % colormap(jet);
    % 
    % 
    % xlabel('$x(t)$ (meters)');
    % ylabel('$y(t)$ (meters)');
    % title('2D Trajectory with Time Gradient');
    % close(v);
 end



function dphi0 = calcDphi0(a, xi0, omega, velocityp, index)
   dphi0 = (velocityp(index))/(xi0-a*(omega(index)*0-(0)));
           
end

function Xdot = diffeq(t, X, g, omega, a, xi0, phi0,index)
    phi = X(1);
    phidot = X(2);
    Xdot = [phidot; ((2*a*omega(index)*phidot)-(a*phidot^2)-(g*sin(phi)))/(xi0 + a*((phi-phi0)- t*omega(index)))];
end

function [value,terminate,direction] = stopfun(t,X, xi0, m, g, phi0, omega, a, index)
    T = m*(g*cos(X(1))+ X(2)^2 * (xi0 - a*(omega(index)*t - (X(1)- phi0))));
    xi = xi0-omega(index)*t*a + (X(1)- phi0)*a;
    value = [xi, T];
    terminate =[1,1];
    direction = [-1, -1];

end