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


close all
clear all
clc

inicond = [2.5556, 0,1.1319,0];
m = [0.40,0.48];
L = [0.197,0.247];
R = [L(1)/2,L(2)/2];
g = 9.81;
I = [0.001839,0.003214];

tspan = [0,30];

Opt = odeset('Stats', 'on', 'Refine', 5, 'RelTol', 1e-5,'AbsTol', 1e-5 ,'Events', [], 'InitialStep', [], 'MaxStep', 0.01);
[t, X] = ode23s(@(t, X) ddiff_eq(t, X, g, L(1), m(1), m(2), I(1), I(2), R(1), R(2)), tspan, inicond, Opt);

           

video_name = sprintf('l.mp4'); 
v = VideoWriter(video_name, 'MPEG-4');

v.FrameRate = 30;  

open(v); 


theta1 = X(:,1); 
theta2 = X(:,3); 

L1 = L(1);
L2 = L(2); 


x1 = L1 * sin(theta1);
y1 = -L1 * cos(theta1);

x2 = x1 + L2 * sin(theta2);
y2 = y1 - L2 * cos(theta2);


dtheta1 = X(:,2);
dtheta2 = X(:,4);


vx1 = L1 * dtheta1 .* cos(theta1);
vy1 = L1 * dtheta1 .* sin(theta1);

vx2 = vx1 + L2 * dtheta2 .* cos(theta2);
vy2 = vy1 + L2 * dtheta2 .* sin(theta2);


velocitySquared1 = sqrt(vx1.^2 + vy1.^2);
velocitySquared2 = sqrt(vx2.^2 + vy2.^2);


interpolation = 100;

t_fine = linspace(t(1), t(end), 1801*interpolation);

[t_unique, ia] = unique(t);
x1_unique = x1(ia);
y1_unique = y1(ia);
x2_unique = x2(ia);
y2_unique = y2(ia);
th1_unique = theta1(ia);
th2_unique = theta2(ia);
dth1_unique = dtheta1(ia);
dth2_unique = dtheta2(ia);

x1_interp = interp1(t_unique, x1_unique, t_fine, 'spline')';
y1_interp = interp1(t_unique, y1_unique, t_fine, 'spline')';
x2_interp = interp1(t_unique, x2_unique, t_fine, 'spline')';
y2_interp = interp1(t_unique, y2_unique, t_fine, 'spline')';


th1_interp = interp1(t_unique, th1_unique, t_fine, 'spline')';
th2_interp = interp1(t_unique, th2_unique, t_fine, 'spline')';
dth1_interp = interp1(t_unique, dth1_unique, t_fine, 'spline')';
dth2_interp = interp1(t_unique, dth2_unique, t_fine, 'spline')';

vs1_unique = velocitySquared1(ia);
vs2_unique = velocitySquared2(ia);

velocitySquared1 = interp1(t_unique, vs1_unique, t_fine, 'spline')';
velocitySquared2 = interp1(t_unique, vs2_unique, t_fine, 'spline')';


velocity_norm1 = (velocitySquared1 - min(velocitySquared1)) / (max(velocitySquared1) - min(velocitySquared1)); 
velocity_norm2 = (velocitySquared2 - min(velocitySquared2)) / (max(velocitySquared2) - min(velocitySquared2)); 


screen_size = get(0, 'ScreenSize');
width = screen_size(3);
height = screen_size(4);

totalfig1 = figure(...
    'Position', [screen_size(3)*0.05, screen_size(4)*0.05, screen_size(3)*0.9, screen_size(4)*0.85], ...
    'Color', 'white', ... 
    'NumberTitle', 'off', ... 
    'Name', 'Part2,3' ...
);

sgtitle('Numerical Results', 'FontSize', 20);

Angles = subplot('Position', [0.1 0.55 0.25 0.35],...
    'Color', 'white');
    
    hold on;
    plot(t_fine, th1_interp, 'Color', [34 / 255, 139 / 255, 34 / 255]); 
    plot(t_fine, th2_interp, 'Color', [0, 0.2, 0.5]); 
    
    grid on;
    
    legend('$\theta_{1}$', '$\theta_{2}$');
    title('$\theta_{1}$ and $\theta_{2}$ in function of time for $t = [0,30]$');
    
    xlabel('time (t)', 'Interpreter', 'latex');
    ylabel('$\theta(rad)$', 'Interpreter', 'latex');
    
    ylim([-30, 70]);  

Acc =  subplot('Position', [0.7 0.55 0.25 0.35],...
    'Color', 'white');

    hold on;
    D = (2*I(2)*L(1)^2*m(2) + 2*I(2)*m(1)*R(1)^2 + L(1)^2*m(2)^2*R(2)^2 + 2*m(1)*m(2)*R(1)^2*R(2)^2 + 2*I(1)*(I(2) + m(2)*R(2)^2) - L(1)^2*m(2)^2*R(2)^2.*cos(2*(X(:,1) - X(:,3))));
    plot(t, (-1*(2*g*m(1)*R(1)*(I(2) + m(2)*R(2)^2).*sin(X(:,1)) + L(1)*m(2)*(g*(2*I(2) + m(2)*R(2)^2).*sin(X(:,1))+ R(2)*(g*m(2)*R(2).*sin(X(:,1) -2*X(:,3)) + 2.*(X(:,4).^2 *(I(2) + m(2)*R(2)^2) + X(:,2).^2*L(1)*m(2)*R(2).*cos(X(:,1) - X(:,3))).*sin(X(:,1) - X(:,3))))))./D, 'Color', [0.9, 0.2, 0]); 
    plot(t, (m(2)*R(2).*(-1*(g*(2*I(1) + L(1)^2*m(2) + 2*m(1)*R(1)^2).*sin(X(:,3))) + L(1)*(g*m(1)*R(1).*sin(X(:,3)) + 2*X(:,2).^2*(I(1) + L(1)^2*m(2) + m(1)*R(1)^2).*sin(X(:,1) - X(:,3))+ X(:,4).^2.*L(1)*m(2)*R(2).*sin(2.*(X(:,1) - X(:,3))) + g*m(1)*R(1).*sin(2*X(:,1) - X(:,3)) + g*L(1)*m(2)*sin(2*X(:,1) - X(:,3)))))./D, 'Color', [0, 0.2, 0.5]); 
    
    grid on;
    
    legend('$\ddot\theta_{1}$', '$\ddot\theta_{2}$');
    title('$\ddot\theta_{1}$ and $\ddot\theta_{2}$ in function of time for $t = [0,30]$');
    
    xlabel('time (t)', 'Interpreter', 'latex');
    ylabel('$\ddot\theta(rad/s_{2})$', 'Interpreter', 'latex');


PhasePortratit1 = subplot('Position', [0.05 0.1 0.4 0.35],...
    'Color', 'white');
    
    hold on;
    plot(X(:,1), X(:,2), 'Color', [0, 0.45, 0.74]); 
 
    
    grid on;
    
    legend('$\dot\theta_1 vs \theta_1$');
    title('Phase portratit 1');
    
    xlabel('$\theta_1(rad)$', 'Interpreter', 'latex');
    ylabel('$\dot\theta_1(rad)$', 'Interpreter', 'latex');
    
PhasePortratit2 = subplot('Position', [0.55 0.1 0.4 0.35],...
    'Color', 'white');
    
    hold on;
 
    plot(X(:,3), X(:,4), 'Color', [0.85, 0.33, 0.10]); 
    
    grid on;
    
    legend('$\dot\theta_2 vs \theta_2$');
    title('Phase portratit 2');
    
    xlabel('$\theta_2(rad)$', 'Interpreter', 'latex');
    ylabel('$\dot\theta_2(rad)$', 'Interpreter', 'latex');
    


Energy = subplot('Position', [0.4 0.55 0.25 0.35],...
    'Color', 'white');

[positions, velocities] = double_pendulum_kinematics(X, L(1), L(2));

[pos_interp, vel_interp] = interpolate_positions(t, t_fine, positions, velocities);


V = m(1)* g * (pos_interp.com1(:,2)) + m(2)* g * (pos_interp.com2(:,2));

T = 0.5*dth1_interp(:,1).^2 .* I(1) + 0.5*dth2_interp(:,1).^2 .* I(2) + 0.5*m(1)*(vel_interp.v_com1(:,1).^2 + vel_interp.v_com1(:,2).^2) + 0.5*m(2)*(vel_interp.v_com2(:,1).^2 + vel_interp.v_com2(:,2).^2);

E = V+T;

    hold on;
    plot(t_fine, V, 'Color', [34 / 255, 139 / 255, 34 / 255]); 
    plot(t_fine, T, 'Color', [0, 0.2, 0.5]); 
    plot(t_fine, E, 'Color', [0.9, 0.2, 0]);
    
    grid on;
    
    legend('$V$', '$T$', '$E$');
    title('$V$, $T$ and $E$ in function of time for $t = [0,30]$');
    
    xlabel('time (t)', 'Interpreter', 'latex');
    ylabel('$E(J)$', 'Interpreter', 'latex');
    ylim([-4, 5]);


    

totalfig2 = figure(...
    'Position', [screen_size(3)*0.05, screen_size(4)*0.05, screen_size(3)*0.9, screen_size(4)*0.825], ...
    'Color', 'white', ... 
    'NumberTitle', 'off', ... 
    'Name', 'Part4' ...
);

sgtitle('Comparison', 'FontSize', 20);

load('tracked_positions.mat');


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


    exp_theta1 = atan2(delta_x1(i), -delta_y1(i)); 
    exp_theta2 = atan2(delta_x2(i), -delta_y2(i)); 


    delta_theta1 = exp_theta1 - prev_theta1;


    if delta_theta1 > pi
        delta_theta1 = delta_theta1 - 2*pi; 
    elseif delta_theta1 < -pi
        delta_theta1 = delta_theta1 + 2*pi;  
    end
    

    delta_theta2 = exp_theta2 - prev_theta2;

  
    if delta_theta2 > pi
        delta_theta2 = delta_theta2 - 2*pi; 
    elseif delta_theta2 < -pi
        delta_theta2 = delta_theta2 + 2*pi;  
    end

    theta1_cont(i,1) = theta1_cont(i-1) + delta_theta1;
    theta2_cont(i,1) = theta2_cont(i-1) + delta_theta2;

    
    prev_theta1 = exp_theta1; 
    prev_theta2 = exp_theta2; 
end


[t_unique, ia] = unique(linspace(t(1), t(end), 1801));
theta1_cont_unique = theta1_cont(ia);
theta2_cont_unique = theta2_cont(ia);
theta1_cont_interp = interp1(t_unique, theta1_cont_unique, t_fine, 'spline')';
theta2_cont_interp = interp1(t_unique, theta2_cont_unique, t_fine, 'spline')';
[t_unique, ia] = unique(t);



DiffAngle1 = subplot('Position', [0.05 0.55 0.4 0.35],...
    'Color', 'white');

    hold on;
    plot(t_fine, th1_interp, 'Color', [0, 0.2, 0.5]); 
    plot(t_fine, theta1_cont_interp, 'Color', [0.9, 0.2, 0]);

    grid on;

    legend('$\theta_{1num}$', '$\theta_{1exp}$');
    title('$\theta_{1exp}$ vs $\theta_{1num}$in function of time for $t = [0,30]$');

    xlabel('time (t)', 'Interpreter', 'latex');
    ylabel('$\theta_1$', 'Interpreter', 'latex');
     ylim([-25, 10]);

DiffAngle2 = subplot('Position', [0.5 0.55 0.4 0.35],...
    'Color', 'white');

    hold on;
    plot(t_fine, th2_interp, 'Color', [0, 0.2, 0.5]); 
    plot(t_fine, theta2_cont_interp, 'Color', [0.9, 0.2, 0]);

    grid on;

    legend('$\theta_{2num}$', '$\theta_{2exp}$');
    title('$\theta_{2exp}$ vs $\theta_{2num}$in function of time for $t = [0,30]$');

    xlabel('time (t)', 'Interpreter', 'latex');
    ylabel('$\theta_2$', 'Interpreter', 'latex');
     ylim([-30, 60]);

DiffPos1 = subplot('Position', [0.05 0.1 0.4 0.35],...
    'Color', 'white');

    hold on;
    plot(x_pivot_1, y_pivot_1, 'Color', [0.1, 0.1, 0.9]); 
    plot(x_pivot_2, y_pivot_2, 'Color', [0.9, 0.1, 0.1]); 
    plot(x_pivot_3, y_pivot_3, 'Color', [0.1, 0.9, 0.1]);

    grid on;

    legend('$Pivot Blue$', '$Pivot Red$', '$Pivot Green$');
    title('Experimental position for $t = [0,30]$');

    xlabel('x(m)', 'Interpreter', 'latex');
    ylabel('$y(m)$', 'Interpreter', 'latex');
    ylim([-0.5, 0.5]);
    xlim([-0.7, 0.7]);


DiffPos2 = subplot('Position', [0.5 0.1 0.4 0.35],...
    'Color', 'white');

    hold on;
    plot(pos_interp.p0(:,1), pos_interp.p0(:,2), 'Color', [0.1, 0.1, 0.9]); 
    plot(pos_interp.p1(:,1), pos_interp.p1(:,2), 'Color', [0.9, 0.1, 0.1]); 
    plot(pos_interp.p2(:,1), pos_interp.p2(:,2), 'Color', [0.1, 0.9, 0.1]);

    grid on;

    legend('$Pivot Blue$', '$Pivot Red$', '$Pivot Green$');
    title('Numerical position for $t = [0,30]$');

    xlabel('x(m)', 'Interpreter', 'latex');
    ylabel('$y(m)$', 'Interpreter', 'latex');
    ylim([-0.5, 0.5]);
    xlim([-0.7, 0.7]);


   
totalfig3 = figure(...
    'Position', [screen_size(3)*0.05, screen_size(4)*0.05, screen_size(3)*0.9, screen_size(4)*0.8], ...
    'Color', 'white', ... 
    'NumberTitle', 'off', ... 
    'Name', 'Poncaire' ...
);

sgtitle('Poncaire Map', 'FontSize', 20);

poncaire_flag = 1;

poncaire_indices = [0,0];

filename_flag = 'poncaire.mat';

if exist(filename_flag, 'file') == 2
    load('poncaire.mat');
end

if poncaire_flag == 1 

    V_const = 0.847781227578096;

    theta1_range = linspace(-pi, pi, 50000); % Rango de -π a π
    theta2_range = linspace(-pi, pi, 5000);
  
    theta1_sol = [];
    theta2_sol = [];

    for angle1 = theta1_range
        for angle2 = theta2_range
        
            V = m(1) * g * ((L(1) / 2) * cos(angle1)) + m(2) * g * ((L(2) / 2) * cos(angle2));
           
            if abs(V - V_const) < 1e-4 
                theta1_sol = [theta1_sol; angle1];
                theta2_sol = [theta2_sol; angle2];
            end
        end
    end

   % that add to the poncaire mat
    n_iter = 5;                       
    
                    
    a = size(poncaire_indices, 1);
    for j = 1:length(t)-1
            
          if X(j,1) < 0 && X(j+1,1) > 0 
                poncaire_indices(a,1) = (X(j,3)+X(j+1,3))/2;
                poncaire_indices(a,2) = (X(j,4)+X(j+1,4))/2;
                a = a+1;
          end
    
    
    end
    
    for i = 1:n_iter
    
        inicond= [theta1_sol(i), 0, theta2_sol(i), 0];
        [t_p, X] = ode23s(@(t_p, X) ddiff_eq(t_p, X, g, L(1), m(1), m(2), I(1), I(2), R(1), R(2)), tspan, inicond, Opt);
       
        for j = 1:length(t_p)-1
            
            if X(j,1) < 0 && X(j+1,1) > 0 
                poncaire_indices(a,1) = (X(j,3)+X(j+1,3))/2;
                poncaire_indices(a,2) = (X(j,4)+X(j+1,4))/2;
                a = a+1;
            end
    
    
        end
        disp(i);
    
    end

    save('poncaire.mat', 'poncaire_indices');
end

    load('poncaire.mat');
    hold on;
 
    scatter(poncaire_indices(:,1), poncaire_indices(:,2), 10, 'MarkerEdgeColor', [0.85, 0.33, 0.10], 'MarkerFaceColor', [0.85, 0.33, 0.10]);
     

    grid on;
    
    legend('$\dot\theta_2 vs \theta_2$');
    
    xlabel('$\theta_2(rad)$', 'Interpreter', 'latex');
    ylabel('$\dot\theta_2(rad)$', 'Interpreter', 'latex');
    ylim([-15, 10]);
    xlim([-5, 5]);


totalfig4 = figure(...
    'Position', [screen_size(3)*0.05, screen_size(4)*0.05, screen_size(3)*0.9, screen_size(4)*0.775], ...
    'Color', 'white', ... 
    'NumberTitle', 'off', ... 
    'Name', 'Video figure' ...
);

sgtitle('Video', 'FontSize', 20);
 
interpolation = 0.05;
% Fine time vector for interpolation
t_fine = linspace(t(1), t(end), length(t)*interpolation);
[t_unique, ia] = unique(t);

x1_interp = interp1(t_unique, x1_unique, t_fine, 'spline')';
y1_interp = interp1(t_unique', y1_unique, t_fine, 'spline')';
x2_interp = interp1(t_unique, x2_unique, t_fine, 'spline')';
y2_interp = interp1(t_unique, y2_unique, t_fine, 'spline')';
th1_interp = interp1(t_unique, th1_unique, t_fine, 'spline')';
th2_interp = interp1(t_unique, th2_unique, t_fine, 'spline')';
dth1_interp = interp1(t_unique, dth1_unique, t_fine, 'spline')';
dth2_interp = interp1(t_unique, dth2_unique, t_fine, 'spline')';


for i = 1:length(t_fine)
    Xlim(i) = 30*i/(length(t)*interpolation);
    xcord2(i,:)= [x1_interp(i)+0.025*cos(th2_interp(i)) x1_interp(i)-0.025*cos(th2_interp(i))  x2_interp(i)-0.025*cos(th2_interp(i))  x2_interp(i)+0.025*cos(th2_interp(i))];
    ycord2(i,:) = [y1_interp(i)+0.025*sin(th2_interp(i))    y1_interp(i)-0.025*sin(th2_interp(i))  y2_interp(i)-0.025*sin(th2_interp(i)) y2_interp(i)+0.025*sin(th2_interp(i))];
    xcord1(i,:)= [0+0.025*cos(th1_interp(i))  0-0.025*cos(th1_interp(i)) x1_interp(i)-0.025*cos(th1_interp(i)) x1_interp(i)+0.025*cos(th1_interp(i))];
    ycord1(i,:) = [0+0.025*sin(th1_interp(i))  0-0.025*sin(th1_interp(i)) y1_interp(i)-0.025*sin(th1_interp(i)) y1_interp(i)+0.025*sin(th1_interp(i))];

    phi = linspace(0, 2*pi, 100);
    r = 0.025;
    x_circle = r * cos(phi);
    y_circle = r * sin(phi);
end

for i = 1:length(t_fine)
    clf;

    hold on;
    subplot('Position', [0.7 0.55 0.25 0.35]);

    hold on;
    plot(t_fine(1:i), th1_interp(1:i), 'Color', [34 / 255, 139 / 255, 34 / 255]); 
    plot(t_fine(1:i), th2_interp(1:i), 'Color', [0, 0.2, 0.5]); 

    grid on;

    legend('$\theta_{1}$', '$\theta_{2}$');
    title('$\theta_{1}$ and $\theta_{2}$ in function of time for $t = [0,30]$');

    xlabel('time (t)', 'Interpreter', 'latex');
    ylabel('$\theta(rad)$', 'Interpreter', 'latex');

    xlim([0,Xlim(i)])

    subplot('Position', [0.7 0.1 0.25 0.35]);

    hold on;
    plot(th1_interp(1:i), dth1_interp(1:i), 'Color', [0, 0.45, 0.74]); 
    plot(th2_interp(1:i), dth2_interp(1:i), 'Color', [0.85, 0.33, 0.10]); 

    grid on;

    legend('$1$', '$2$');
    title('Phase portratit');

    xlabel('$\theta(rad)$', 'Interpreter', 'latex');
    ylabel('$\dot\theta(rad)$', 'Interpreter', 'latex');



    subplot('Position', [0.05 0.1 0.6 0.8]);

    grid on;
    hold on;

    plot(0.7, 0.7, 'bo', 'MarkerSize', 0.08, 'MarkerFaceColor', 'r'); % Masa 1
    plot(-0.7, -0.7, 'bo', 'MarkerSize', 0.08, 'MarkerFaceColor', 'r'); % Masa 1

    xlabel('$x(t)$ (meters)', 'Interpreter', 'latex');
    ylabel('$z(t)$ (meters)', 'Interpreter', 'latex');
    title(sprintf('3D Video with velocity color gradient'), 'Interpreter', 'latex');

    axis equal;

    fill(xcord2(i,:), ycord2(i,:), 'black'); % 'r' indica color rojo

    fill(xcord1(i,:), ycord1(i,:), 'black'); % 'r' indica color rojo



    fill(x_circle, y_circle, 'black', 'EdgeColor', [0.2 0.2 0.6], 'LineWidth', 1);
    fill(x1_interp(i)+x_circle, y1_interp(i)+y_circle, 'black', 'EdgeColor', [0.2 0.2 0.6], 'LineWidth', 1);
    fill(x2_interp(i)+x_circle, y2_interp(i)+y_circle, 'black', 'EdgeColor', [0.2 0.2 0.6], 'LineWidth', 1);

    plot(0, 0, 'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b'); % Masa 1
    plot(x1_interp(i), y1_interp(i), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r'); % Masa 1
    plot(x2_interp(i), y2_interp(i), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g'); % Masa 2

    elapsed_time = t_fine(i); 
    text(0, 0.6, sprintf('%.1f s', elapsed_time), 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');

    frame = getframe(gcf);
    %writeVideo(v, frame);

    if i > 1
        pause(t_fine(i) - t_fine(i - 1)); 
    end
end
    for i = 1:200

        frame = getframe(gcf);


        %writeVideo(v, frame);
    end
    % colormap(jet);
    close(v);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


function Thetadot = ddiff_eq(t, X, g, L1, m1, m2, I1, I2, R1, R2)

    D = (2*I2*L1^2*m2 + 2*I2*m1*R1^2 + L1^2*m2^2*R2^2 + 2*m1*m2*R1^2*R2^2 + 2*I1*(I2 + m2*R2^2) - L1^2*m2^2*R2^2*cos(2*(X(1) - X(3))));

    Thetadot(1,:) = X(2);
   
    Thetadot(2,:) = (-1*(2*g*m1*R1*(I2 + m2*R2^2)*sin(X(1)) + L1*m2*(g*(2*I2 + m2*R2^2)*sin(X(1))+ R2*(g*m2*R2*sin(X(1) -2*X(3)) + 2*(X(4)^2 *(I2 + m2*R2^2) + X(2)^2*L1*m2*R2*cos(X(1) - X(3)))*sin(X(1) - X(3))))))/D;
    
    Thetadot(3,:) = X(4);

    Thetadot(4,:) = (m2*R2*(-1*(g*(2*I1 + L1^2*m2 + 2*m1*R1^2)*sin(X(3))) + L1*(g*m1*R1*sin(X(3)) + 2*X(2)^2*(I1 + L1^2*m2 + m1*R1^2)*sin(X(1) - X(3))+ X(4)^2*L1*m2*R2*sin(2*(X(1) - X(3))) + g*m1*R1*sin(2*X(1) - X(3)) + g*L1*m2*sin(2*X(1) - X(3)))))/D;
          
end


function [positions, velocities] = double_pendulum_kinematics(X, L1, L2)

    theta1 = X(:,1);
    theta1_dot = X(:,2);
    theta2 = X(:,3);
    theta2_dot = X(:,4);

    x2 = L1 * sin(theta1);
    y2 = -L1 * cos(theta1);

    x3 = x2 + L2 * sin(theta2);
    y3 = y2 - L2 * cos(theta2);

    x_com2 = (L1 / 2) * sin(theta1);
    y_com2 = -(L1 / 2) * cos(theta1);

    x_com3 = x2 + (L2 / 2) * sin(theta2);
    y_com3 = y2 - (L2 / 2) * cos(theta2);

    v_com2_x = (L1 / 2) * theta1_dot .* cos(theta1);
    v_com2_y = (L1 / 2) * theta1_dot .* sin(theta1);

    v2_x = L1 * theta1_dot .* cos(theta1);
    v2_y = L1 * theta1_dot .* sin(theta1);

    v_com3_x = v2_x + (L2 / 2) .* theta2_dot .* cos(theta2);
    v_com3_y = v2_y + (L2 / 2) .* theta2_dot .* sin(theta2);

    x1 = x2*0;
    y1 = x2*0;

 
    positions = struct('p0', [x1, y1], 'p1', [x2, y2], 'p2', [x3, y3], ...
                       'com1', [x_com2, y_com2], 'com2', [x_com3, y_com3]);
    velocities = struct('v_com1', [v_com2_x, v_com2_y], ...
                        'v_com2', [v_com3_x, v_com3_y]);
end

function [pos_interp, vel_interp] = interpolate_positions(t, t_fine, positions, velocities)

    [t_unique, ia] = unique(t);
   
    p01 = positions.p0(:,1)';
    p02 = positions.p0(:,2)';

    p11 = positions.p1(:,1)';
    p12 = positions.p1(:,2)';
    
    p21 = positions.p2(:,1)';
    p22 = positions.p2(:,2)';
    
    com1_1 = positions.com1(:,1)';
    com1_2 = positions.com1(:,2)';
    
    com2_1 = positions.com2(:,1)';
    com2_2 = positions.com2(:,2)';
    
    v_com1_1 = velocities.v_com1(:,1)';
    v_com1_2 = velocities.v_com1(:,2)';
    
    v_com2_1 = velocities.v_com2(:,1)';
    v_com2_2 = velocities.v_com2(:,2)';

    p01_unique = p01(ia);
    p02_unique = p02(ia);
    
    p11_unique = p11(ia);
    p12_unique = p12(ia);
    
    p21_unique = p21(ia);
    p22_unique = p22(ia);
    
    com1_1_unique = com1_1(ia);
    com1_2_unique = com1_2(ia);
    
    com2_1_unique = com2_1(ia);
    com2_2_unique = com2_2(ia);
    
    v_com1_1_unique = v_com1_1(ia);
    v_com1_2_unique = v_com1_2(ia);
    
    v_com2_1_unique = v_com2_1(ia);
    v_com2_2_unique = v_com2_2(ia);

  
    % Perform interpolation
    p01_interp = interp1(t_unique', p01_unique', t_fine, 'spline');
    p02_interp = interp1(t_unique, p02_unique', t_fine, 'spline');
    
    p11_interp = interp1(t_unique', p11_unique', t_fine, 'spline');
    p12_interp = interp1(t_unique', p12_unique', t_fine, 'spline');
    
    p21_interp = interp1(t_unique', p21_unique', t_fine, 'spline');
    p22_interp = interp1(t_unique', p22_unique', t_fine, 'spline');
    
    com1_1_interp = interp1(t_unique', com1_1_unique', t_fine, 'spline');
    com1_2_interp = interp1(t_unique', com1_2_unique', t_fine, 'spline');
    
    com2_1_interp = interp1(t_unique', com2_1_unique', t_fine, 'spline');
    com2_2_interp = interp1(t_unique', com2_2_unique', t_fine, 'spline');
    
    v_com1_1_interp = interp1(t_unique', v_com1_1_unique', t_fine, 'spline');
    v_com1_2_interp = interp1(t_unique', v_com1_2_unique', t_fine, 'spline');
    
    v_com2_1_interp = interp1(t_unique', v_com2_1_unique', t_fine, 'spline');
    v_com2_2_interp = interp1(t_unique', v_com2_2_unique', t_fine, 'spline');
    
    % Collect the interpolated values in a structure
    
    pos_interp = struct('p0', [p01_interp; p02_interp]', 'p1', [p11_interp; p12_interp]', 'p2', [p21_interp; p22_interp]', ...
                       'com1', [com1_1_interp; com1_2_interp]', 'com2', [com2_1_interp; com2_2_interp]');
    vel_interp = struct('v_com1', [v_com1_1_interp; v_com1_2_interp]', ...
                        'v_com2', [v_com2_1_interp; v_com2_2_interp]');
end
