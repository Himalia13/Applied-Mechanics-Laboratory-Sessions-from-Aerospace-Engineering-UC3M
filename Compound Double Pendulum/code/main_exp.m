% % Numerical Solutions

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

close all; clear; clc;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                           PLAY WITH THIS:
flag.get_RF = 1; 

Lp1 = 0.197; 
Lp2 = 0.247; 
m = [0.40,0.48];
g = 9.81;
I = [0.001839,0.003214];

load('tracked_positions.mat');


theta1_cont = zeros(1, length(x_pivot_1));  
theta2_cont = zeros(1, length(x_pivot_1));  
time = linspace(0, 30, length(x_pivot_1));      


delta_x1 = x_pivot_2 - x_pivot_1;
delta_y1 = y_pivot_2 - y_pivot_1;
delta_x2 = x_pivot_3 - x_pivot_2;
delta_y2 = y_pivot_3 - y_pivot_2;

theta1_cont(1) = atan2(delta_x1(1), -delta_y1(1)); 
theta2_cont(1) = atan2(delta_x2(1), -delta_y2(1));

fprintf('Initial angles:\n');
fprintf('Theta1: %.4f, Theta2: %.4f\n', theta1_cont(1), theta2_cont(1));

prev_theta1 = theta1_cont(1);  
prev_theta2 = theta2_cont(1);  

for i = 2:length(x_pivot_1)

    delta_x1(i) = x_pivot_2(i) - x_pivot_1(i);
    delta_y1(i) = y_pivot_2(i) - y_pivot_1(i);

  
    delta_x2(i) = x_pivot_3(i) - x_pivot_2(i);
    delta_y2(i) = y_pivot_3(i) - y_pivot_2(i);


    theta1 = atan2(delta_x1(i), -delta_y1(i));
    theta2 = atan2(delta_x2(i), -delta_y2(i));  

    delta_theta1 = theta1 - prev_theta1;

 
    if delta_theta1 > pi
        delta_theta1 = delta_theta1 - 2*pi; 
    elseif delta_theta1 < -pi
        delta_theta1 = delta_theta1 + 2*pi; 
    end
    
    delta_theta2 = theta2 - prev_theta2;

    if delta_theta2 > pi
        delta_theta2 = delta_theta2 - 2*pi;  
    elseif delta_theta2 < -pi
        delta_theta2 = delta_theta2 + 2*pi;  
    end

  
    theta1_cont(i) = theta1_cont(i-1) + delta_theta1;
    theta2_cont(i) = theta2_cont(i-1) + delta_theta2;

  
    prev_theta1 = theta1;  
    prev_theta2 = theta2; 
end

screen_size = get(0, 'ScreenSize');
width = screen_size(3);
height = screen_size(4);

totalfig = figure(...
    'Position', [screen_size(3)*0.05, screen_size(4)*0.05, screen_size(3)*0.9, screen_size(4)*0.85], ...
    'Color', 'white', ... 
    'NumberTitle', 'off', ... 
    'Name', 'Experimental Part' ...
);

sgtitle('Experimental Part', 'FontSize', 20);


Angles = subplot('Position', [0.1 0.1 0.25 0.35],...
    'Color', 'white');
    
    hold on;
    plot(time, theta1_cont, 'Color', [34 / 255, 139 / 255, 34 / 255]); 
    plot(time, theta2_cont, 'Color', [0, 0.2, 0.5]); 
    
    grid on;
    
    legend('$\theta_{1}$', '$\theta_{2}$');
    title('$\theta_{1}$ and $\theta_{2}$ in function of time for $t = [0,30]$');
    
    xlabel('time (t)', 'Interpreter', 'latex');
    ylabel('$\theta(rad)$', 'Interpreter', 'latex');

Energy = subplot('Position', [0.4 0.1 0.25 0.35],...
    'Color', 'white');

dt = 1/60;
omega1(1) = 0;
omega2(1) = 0;
N=length(x_pivot_1);
  

    for i = 1:N-1
        rx1 = x_pivot_2(i) - x_pivot_1(i); 
        ry1 = y_pivot_2(i) - y_pivot_1(i);
        rx_next1 = x_pivot_2(i+1)- x_pivot_1(i+1);
        ry_next1 = y_pivot_2(i+1) - y_pivot_1(i+1);

    
        dot_product1 = rx1 * rx_next1 + ry1 * ry_next1;
        magnitude_r1 = sqrt(rx1^2 + ry1^2);
        magnitude_r1next = sqrt(rx_next1^2 + ry_next1^2);
        
   
        cos_theta1 = dot_product1 / (magnitude_r1 * magnitude_r1next);
        delta_theta1 = acos(cos_theta1);
        
  
        cross_product1 = rx1 * ry_next1 - ry1 * rx_next1;
        if cross_product1 < 0
            delta_theta1 = -delta_theta1; 
        end
        
        omega1(i+1) = delta_theta1 / dt; 
    end
    omega1 = omega1';
    for i = 1:N-1
        rx2 = x_pivot_3(i) - x_pivot_2(i); 
        ry2 = y_pivot_3(i) - y_pivot_2(i);
        rx_next2 = x_pivot_3(i+1)- x_pivot_2(i+1);
        ry_next2 = y_pivot_3(i+1) - y_pivot_2(i+1);

        dot_product2 = rx2 * rx_next2 + ry2 * ry_next2;
        magnitude_r2 = sqrt(rx2^2 + ry2^2);
        magnitude_r2next = sqrt(rx_next2^2 + ry_next2^2);
        
    
        cos_theta2 = dot_product2 / (magnitude_r2 * magnitude_r2next);
        delta_theta2 = acos(cos_theta2);
        
    
        cross_product2 = rx2 * ry_next2 - ry1 * rx_next2;
        if cross_product2 < 0
            delta_theta2 = -delta_theta2; % Sentido horario
        end
        
 
        omega2(i+1) = delta_theta2 / dt; 
    end
    omega2 = omega2';


    vG1(1)= 0;
    

   
        Gx1 = (x_pivot_1 + x_pivot_2) / 2; % (x1 + x2) / 2
        Gy1 = (y_pivot_1 + y_pivot_2) / 2; % (y1 + y2) / 2
   
    
 
    for i = 1:N-1
        dGx1 = (Gx1(i+1) - Gx1(i)) / dt; 
        dGy1 = (Gy1(i+1) - Gy1(i)) / dt; 
        vG1(i+1) = sqrt(dGx1^2 + dGy1^2);  
    end
        vG1 = vG1';
   
      vG2(1)= 0;
    

 
        Gx2 = (x_pivot_2 + x_pivot_3) / 2;
        Gy2 = (y_pivot_2 + y_pivot_3) / 2;
   
    
  
    for i = 1:N-1
        dGx2 = (Gx2(i+1) - Gx2(i)) / dt; 
        dGy2 = (Gy2(i+1) - Gy2(i)) / dt; 
        vG2(i+1) = sqrt(dGx2^2 + dGy2^2);  
    end
    vG2 = vG2';

V = m(1)* g * (Gy1) + m(2)* g * ((Gy2));

T = 0.5*omega1.^2 .* I(1) + 0.5*omega2.^2 .* I(2) + 0.5*m(1)*vG1.^2 + 0.5*m(2)*vG2.^2;
E = V+T;

    hold on;
    plot(time, V, 'Color', [34 / 255, 139 / 255, 34 / 255]); 
    plot(time, T, 'Color', [0, 0.2, 0.5]); 
    plot(time, E, 'Color', [0.9, 0.2, 0]);
    
    grid on;
    
    legend('$V$', '$T$', '$E$');
    title('$V$, $T$ and $E$ in function of time for $t = [0,30]$');
    
    xlabel('time (t)', 'Interpreter', 'latex');
    ylabel('$E(J)$', 'Interpreter', 'latex');

CM_Velocity = subplot('Position', [0.7 0.1 0.25 0.35],...
    'Color', 'white');
    
    hold on;
    plot(time, vG1, 'Color', [34 / 255, 139 / 255, 34 / 255]); 
    plot(time, vG2, 'Color', [0, 0.2, 0.5]); 
    
    grid on;
    
    legend('$\vec{v}_{G1}$', '$\vec{v}_{G2}$');
    title('$\vec{v}_{G1}$ and $\vec{v}_{G2}$ in function of time for $t = [0,30]$');
    
    xlabel('time (t)', 'Interpreter', 'latex');
    ylabel('$\vec{v}_{G}(m/s)$', 'Interpreter', 'latex');


Pivot2 = subplot('Position', [0.1 0.55 0.25 0.35],...
    'Color', 'white');

    plot(x_pivot_2, y_pivot_2, 'r', 'LineWidth', 0.5);
    axis equal;
    axis([referenceFrame.xAxis(1), referenceFrame.xAxis(end), ...
        referenceFrame.yAxis(end), referenceFrame.yAxis(1)]);
    hold on;
   
    hold off;
    title('$X_{red}$ and $Y_{red}$ for $t = [0,30]$');
    xlabel('X (m)');
    ylabel('Y (m)');
    grid on; grid minor;


Pivot3 = subplot('Position', [0.4 0.55 0.25 0.35],...
    'Color', 'white');
    plot(x_pivot_3, y_pivot_3, 'g', 'LineWidth', 0.5);
    axis equal;
    axis([referenceFrame.xAxis(1), referenceFrame.xAxis(end), ...
        referenceFrame.yAxis(end), referenceFrame.yAxis(1)]);
    hold on;
  
    title('$X_{green}$ and $Y_{green}$ for $t = [0,30]$');
    xlabel('X (m)');
    ylabel('Y (m)');
    grid on; grid minor;


Total = subplot('Position', [0.7 0.55 0.25 0.35],...
    'Color', 'white');
    plot(x_pivot_3-x_pivot_2, y_pivot_3-y_pivot_2, 'k', 'LineWidth', 0.5);
    axis equal;
    axis([-1, 1, ...
        -1, 1]*Lp2*1.2);
    hold on;
    
    title('$X_{green}-X_{red}$ and $Y_{green}-Y_{red}$ for $t = [0,30]$');
    xlabel('$X_{green}-X_{red}$ (m)');
    ylabel('$Y_{green}-Y_{red}$ (m)');
    grid on; grid minor;


