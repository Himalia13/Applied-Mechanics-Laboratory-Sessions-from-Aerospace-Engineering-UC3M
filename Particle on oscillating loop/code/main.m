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

A = pi;
g = 9.81; % m/s`2
m = 1; % kg
R = 1; % m

omega = 0.1;
phi0 = 0;
dphi0 = 0;

    
tspan = [0 3*2* pi /omega];
X0 = [phi0; dphi0];

%solving the ODE with ode45
Opt = odeset('Stats', 'on', 'Refine', 20, 'RelTol', 1e-8,'AbsTol', 1e-10 , 'MaxStep', 0.01);%'Refine', 10
[t, X] = ode45(@(t, X) diffeq(t, X, g, omega, A, R), tspan, X0, Opt);

omegaT = A * omega * sin(omega*t);
theta = -A * cos(omega*t) + A;



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



%Kinematics

x = R.*(1 + cos(X(:,1))) .* cos(theta);
y = R.*(1 + cos(X(:,1))) .* sin(theta);
z = R.*(sin(X(:,1)));

velocity = [R.*X(:,2).*(-sin(X(:,1)).*cos(theta)) + omegaT.*(1+cos(X(:,1))).*(-sin(theta)), R.*X(:,2).*cos(X(:,1)), R.*X(:,2).*(-sin(X(:,1)).*sin(theta))+omegaT.*(1+cos(X(:,1))).*(-sin(theta))];

velocitySquared = sqrt(velocity(:,1).^2 + velocity(:,3).^2 + velocity(:,2).^2);

%ENERGY

%Potential Energy V = mg(y-yo)
posX0 = 0;

V = m*g*(z - posX0)+ 0.5 * m.* (omegaT).^2 .* (x.^2 + y.^2);

%Kinetic Energy T = 0.5m*v*2
T = 0.5*m*(velocitySquared.^2);
%Mechanic Energy E = V+T
Energy=V+T;

%Plot the Energy as a function of time.

fig1 = figure(1);
set(fig1, 'Position', [500, 500, 560, 420]);
hold on


plot(t, Energy, 'Color', [139 / 255, 0 / 255, 0 / 255]); 
plot(t, V, 'Color', [34 / 255, 139 / 255, 34 / 255]); 
plot(t, T, 'Color', [0, 0.2, 0.5]); 

grid on;

legend('$E(J)$', 'Potential', 'Kinetic');
title('Energy in function of time for $t = [0,T]$', 'Interpreter', 'latex');

xlabel('time (t)', 'Interpreter', 'latex');
ylabel('$E(J)$', 'Interpreter', 'latex');

ylim([-20, 20]);
xlim([0, 2* pi /omega]);

% print(1, '-dpng', '-r600');

fig1_2 = figure(12);
set(fig1_2, 'Position', [1060, 500, 560, 420]);
hold on

plot(t, Energy, 'Color', [139 / 255, 0 / 255, 0 / 255]); 
plot(t, V, 'Color', [34 / 255, 139 / 255, 34 / 255]); 
plot(t, T, 'Color', [0, 0.2, 0.5]); 

grid on;

legend('$E(J)$', 'Potential', 'Kinetic');
title('Energy in function of time for $t = [0,3T]$', 'Interpreter', 'latex');

xlabel('time (t)', 'Interpreter', 'latex');
ylabel('$E(J)$', 'Interpreter', 'latex');

ylim([-20, 20]);
xlim([0, 3*2* pi /omega]);

% print(12, '-dpng', '-r600');



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



%Reaction Forces

N1 = m * g * sin(X(:, 1)) - m * A^2 * R * omega^2 * (1 + cos(X(:, 1))) .* sin(omega * t) .^ 2 .* cos(X(:, 1)) - m * R * X(:, 2).^2;

N2 = 2 * m * A * R * X(:, 2) .* sin(X(:, 1)) .* sin(omega * t) - m * A * R * omega^2 * (1 + cos(X(:, 1))) .* cos(omega * t);

pause(0.5)
fig2 = figure(2);
set(fig2, 'Position', [500, 500, 560, 420]);
hold on


plot(t, N1, 'Color', [0, 0.45, 0.74]); 
plot(t, N2, 'Color', [0.85, 0.33, 0.10]); 

grid on;

legend('$N_{1}$', '$N_{2}$', 'Location', 'north');
title('Reaction forces $\vec{N} = N_{1} \vec{e_{r}} \quad N_{2}\vec{k}$ for $t = [0,T]$', 'Interpreter', 'latex');

xlabel('time (t)', 'Interpreter', 'latex');
ylabel('$N_{1,2}(N)$', 'Interpreter', 'latex');

ylim([-30, 30]);
xlim([0, 2* pi /omega]);

% print(2, '-dpng', '-r600');

fig2_2 = figure(22);
set(fig2_2, 'Position', [1060, 500, 560, 420]);
hold on

plot(t, N1, 'Color', [0, 0.45, 0.74]); 
plot(t, N2, 'Color', [0.85, 0.33, 0.10]);

grid on;

legend('$N_{1}$', '$N_{2}$', 'Location', 'north');
title('Reaction forces $\vec{N} = N_{1} \vec{e_{r}} \quad N_{2}\vec{k}$ for $t = [0,3T]$ ', 'Interpreter', 'latex');

xlabel('time (t)', 'Interpreter', 'latex');
ylabel('$N_{1,2}(N)$', 'Interpreter', 'latex');

ylim([-50, 50]);
xlim([0, 3*2* pi /omega]);

% print(22, '-dpng', '-r600');



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



F_Centrip = m * A^2 * R * omega^2 * (1 + cos(X(:,1))) .* (sin(omega * t)).^2;

F_Corio = -2 * m * A * R * omega .* X(:,2) .* sin(X(:,1)) .* sin(omega * t);

F_Euler = m * A * R * omega^2 .* (1 + cos(X(:,1))) .* cos(omega * t);

pause(0.5)
fig3 = figure(3);
set(fig3, 'Position', [500, 500, 560, 420]);
hold on


plot(t, F_Centrip, 'Color', [139 / 255, 0 / 255, 0 / 255]); 
plot(t, F_Corio, 'Color', [34 / 255, 139 / 255, 34 / 255]); 
plot(t, F_Euler, 'Color', [0, 0.2, 0.5]); 

grid on;

legend('$\vec{F}_{centripelal}$', '$\vec{F}_{coriolis}$', '$\vec{F}_{euler}$', 'Location', 'east', 'Location', 'north');
title('Fictitious Forces  for $t = [0,T]$', 'Interpreter', 'latex');

xlabel('time (t)', 'Interpreter', 'latex');
ylabel('$F(N)$', 'Interpreter', 'latex');

ylim([-3, 3]);
xlim([0, 2* pi /omega]);

% print(3, '-dpng', '-r600');

fig3_2 = figure(32);
set(fig3_2, 'Position', [1060, 500, 560, 420]);
hold on

plot(t, F_Centrip, 'Color', [139 / 255, 0 / 255, 0 / 255]); 
plot(t, F_Corio, 'Color', [34 / 255, 139 / 255, 34 / 255]); 
plot(t, F_Euler, 'Color', [0, 0.2, 0.5]); 

grid on;

legend('$\vec{F}_{centripelal}$', '$\vec{F}_{coriolis}$', '$\vec{F}_{euler}$', 'Location', 'east', 'Location', 'north');
title('Fictitious Forces  for $t = [0,3T]$', 'Interpreter', 'latex');

xlabel('time (t)', 'Interpreter', 'latex');
ylabel('$F(N)$', 'Interpreter', 'latex');

ylim([-5, 5]);
xlim([0, 3*2* pi /omega]);

% print(32, '-dpng', '-r600');



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



% Plotting the trajectory with color gradient
pause(0.5)
t_fine = linspace(min(t), max(t), 50000); % More points for smooth curve

velocitySquared = sqrt(velocity(:,1).^2 + velocity(:,3).^2 + velocity(:,2).^2);

[t_unique, ia] = unique(t);
x_unique = x(ia);
y_unique = y(ia);
z_unique = z(ia);
velocitySquared = velocitySquared(ia);
x_interp = interp1(t_unique, x_unique, t_fine, 'spline');
y_interp = interp1(t_unique, y_unique, t_fine, 'spline');
z_interp = interp1(t_unique, z_unique, t_fine, 'spline');
velocitySquared = interp1(t_unique, velocitySquared, t_fine, 'spline');


velocity_norm = (velocitySquared - min(velocitySquared)) / (max(velocitySquared) - min(velocitySquared)); 

cmap = colormap(jet); 
num_colors = size(cmap, 1);

color_indices = round(velocity_norm * (num_colors - 1)) + 1;
colorg1 = cmap(color_indices, :); 


pause(0.5)
fig4 = figure(4);
set(fig4, 'Position', [660, 420, 800, 500]);
hold on;

scatter3(x_interp, y_interp, z_interp, 1, colorg1, 'filled');  
colormap(jet);  

xlabel('$x(t)$ (meters)', 'Interpreter', 'latex');
ylabel('$z(t)$ (meters)', 'Interpreter', 'latex');
zlabel('$y(t)$ (meters)', 'Interpreter', 'latex');
title('3D Trajectory with Velocity Gradient', 'Interpreter', 'latex');

grid on;
axis equal;
view(3);


cb = colorbar;
cb.Label.String = 'Velocity (m/s)'; 
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 14;

cb.Ticks = linspace(min(velocitySquared), max(velocitySquared), 5); 
cb.TickLabels = arrayfun(@(v) sprintf('%.1f', v), cb.Ticks, 'UniformOutput', false);

clim([min(velocitySquared) max(velocitySquared)]); 

% print(4, '-dpng', '-r600'); 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% % % % % Uncomment this section if you want to save the videos.
% % 
% % video_name = sprintf('trajectory_video.mp4'); 
% % v = VideoWriter(video_name, 'MPEG-4');
% % 
% % v.FrameRate = 60;  % frame rate
% % 
% % open(v); 
% % 
% % % Fine time vector for interpolation
% % t_fine = linspace(min(t), max(t), 60*max(t)); % More points for smooth curve
% % 
% % theta_rotation = theta;
% % 
% % 
% % velocitySquared = sqrt(velocity(:,1).^2 + velocity(:,3).^2 + velocity(:,2).^2);
% % 
% % [t_unique, ia] = unique(t);
% % x_unique = x(ia);
% % y_unique = y(ia);
% % z_unique = z(ia);
% % theta_rotation = theta_rotation(ia);
% % x_interp = interp1(t_unique, x_unique, t_fine, 'spline');
% % y_interp = interp1(t_unique, y_unique, t_fine, 'spline');
% % z_interp = interp1(t_unique, z_unique, t_fine, 'spline');
% % theta_rotation = interp1(t_unique, theta_rotation, t_fine, 'spline');
% % [t_unique, ia] = unique(t);
% % X_unique = X(ia, 1);
% % velocitySquared = interp1(t_unique, velocitySquared, t_fine, 'spline');
% % 
% % 
% % velocity_norm = (velocitySquared - min(velocitySquared)) / (max(velocitySquared) - min(velocitySquared)); 
% % 
% %
% % cmap = colormap(jet); 
% % num_colors = size(cmap, 1); 
% % 
% % 
% % color_indices = round(velocity_norm * (num_colors - 1)) + 1; 
% % colorg1 = cmap(color_indices, :); 
% % 
% % 
% % total_time = max(t);
% % 
% % 
% % center = [1, 0, 0]; 
% % 
% % 
% % 
% % figure_handle = figure(5); 
% % set(figure_handle, 'Position', [0, 0, 1920, 1080]);  
% % grid on;
% % hold on;
% % 
% % 
% % for i = 1:length(t_fine)
% %     clf; 
% %     view(3);
% %     grid on;
% %     hold on;
% % 
% %     scatter3(2, 2, 0, 1,'b', 'filled');
% %     hold on;
% % 
% %      scatter3(-2, -2, -1, 1, 'b', 'filled');
% %     hold on;
% % 
% %     scatter3(x_interp(1:i), y_interp(1:i), z_interp(1:i), 5, colorg1(1:i, :), 'filled');%cmap_jet(1:i, :)
% %     hold on;
% % 
% %     scatter3(x_interp(i), y_interp(i), z_interp(i), 60, colorg1(i, :), 'filled');%cmap_jet(i, :)
% % 
% %     xlabel('$x(t)$ (meters)', 'Interpreter', 'latex');
% %     ylabel('$y(t)$ (meters)', 'Interpreter', 'latex');
% %     zlabel('$z(t)$ (meters)', 'Interpreter', 'latex');
% %     title(sprintf('Case %d: 3D Video', 1), 'Interpreter', 'latex');
% % 
% %     axis equal;
% % 
% %     plot3([0 0], [0 0], [1 -1], 'black', 'LineWidth', 2); % Eje X en rojo
% % 
% %     theta_circle = linspace(0, 2*pi, 100); 
% % 
% %     x_circle = center(1) + R * cos(theta_circle);
% %     z_circle = center(3) + R * sin(theta_circle);
% %     y_circle = zeros(size(x_circle));  
% % 
% %     Rot = [cos(theta_rotation(i)), -sin(theta_rotation(i)), 0; sin(theta_rotation(i)),  cos(theta_rotation(i)), 0; 0, 0, 1];
% % 
% %     coords_rotated = Rot * [x_circle; y_circle; z_circle];
% %     hcircle = plot3(coords_rotated(1, :), coords_rotated(2, :),coords_rotated(3, :), 'black', 'LineWidth', .5);
% % 
% %     elapsed_time = t_fine(i); 
% %     text(1.1, 1.5, sprintf('%.1f s', elapsed_time), 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');
% % 
% %     frame = getframe(gcf);
% %     writeVideo(v, frame);
% % 
% %     if i > 1
% %         pause(t_fine(i) - t_fine(i - 1)); 
% %     end
% % end
% %     for i = 1:200
% % 
% %         frame = getframe(gcf);
% % 
% % 
% %         writeVideo(v, frame);
% %     end
% %     % colormap(jet);
% %     close(v);
% % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


pause(0.5)
A = 0;
phi0_values = [0, pi/4, pi/2];  
dphi0_values = [0, 1.5, 3];      

fig6 = figure(6);
set(fig6, 'Position', [360, 420, 700, 500]);
colors = [
    0, 0, 1;     
    0, 0.5, 0;   
    1, 0, 0;   
    0, 1, 1;   
    0.5, 0, 0.5;
    1, 0.5, 0;   
    0.5, 0.5, 0; 
    0, 0, 0;     
    0.5, 0.5, 1  
];

h = zeros(length(phi0_values) * length(dphi0_values), 1); 
legendStrings = cell(length(phi0_values) * length(dphi0_values), 1);


for i = 1:length(phi0_values)
    for j = 1:length(dphi0_values)
        phi0 = phi0_values(i);
        dphi0 = dphi0_values(j);

        tspan = [0 3*2* pi / omega];
        X0 = [phi0; dphi0];
        
        index = (i - 1) * 3 + j;

        switch index 
            case 1
                color = colors(1, :); 
            case 2
                color = colors(2, :); 
            case 3
                color = colors(3, :); 
            case 4
                color = colors(4, :); 
            case 5
                color = colors(5, :);
            case 6
                color = colors(6, :); 
            case 7
                color = colors(7, :);
            case 8
                color = colors(8, :); 
            case 9
                color = colors(9, :); 
            otherwise
                color = [0, 0, 0];
        end
        
        Opt = odeset('Stats', 'on', 'Refine', 30, 'RelTol', 1e-12, 'AbsTol', 1e-14, 'MaxStep', 1e-2);
        if ~isequal(X0, [pi/2;0])
            [t, X] = ode45(@(t, X) diffeq(t, X, g, omega, A, R), tspan, X0, Opt);
            h(index) = plot(X(:, 1), X(:, 2), 'k'); 
            set(h(index), 'Color', color);  
        else
              t = 0;
              X = [X0(1,1),X0(2,1)]; 
              h(index) = plot(X(1,1), X(1, 2), 'o'); 
              set(h(index), 'Color', color);  
            
             
        end
      
        

        fprintf('\nPlease wait, this could take a bit more time.\n')

        hold on
        
        title('$\dot\phi$ in function of $\phi$.', 'Interpreter', 'latex');
        xlabel('$\phi$ (rad)');
        
        
       
        ylabel('$\dot\phi$ (rad)');
        
        legendStrings{index} = sprintf('$\\phi_{0} = %.2f$, $\\dot{\\phi}_{0} = %.2f$', phi0, dphi0);
        plot(X(1,1), X(1,2), 'x', 'Color', color, 'HandleVisibility', 'off');
        
        grid on;
      
    end
end

ylim([-8, 8]);
xlim([-6, 8]);

legend(h, legendStrings, 'Interpreter', 'latex', 'Location', 'best');

% print(6, '-dpng', '-r600');



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


pause(0.5)
A=0.1;
R = 0.01;

omega_values = [sqrt(g / R), (1/2) * sqrt(g / R), (1/4) * sqrt(g / R)];

delta0 = 0;
ddelta0 = 0;

phi0 = -pi/2;
dphi0 = 0;

X0 = [phi0; dphi0];
Y0 = [delta0; ddelta0];

legendStrings2 = cell(3, 1);

fig7 = figure(7);
set(fig7, 'Position', [1060, 420, 700, 500]);
for i = 1:3
    tspan = linspace(0, 3*2* pi /omega_values(i),1000);

    Opt = odeset('Stats','on', 'MaxStep', 1e-4, 'Refine', 30, 'RelTol', 1e-10,'AbsTol', 1e-12);
    
    hold on;
    [t_x, X] = ode45(@(t, X) diffeq(t, X, g, omega_values(i), A, R), tspan, X0, Opt);
    phi  = X(:,1);
  
    [t, Y] = ode45(@(t, Y) diffeq1(t, Y, g, omega_values(i), A, R), tspan, Y0, Opt);
    delta = Y(:,1);

    relative_diff = (phi - delta + (pi/2)) .* 100 ./ phi;
    

    h(i) = plot(t, relative_diff, 'k'); 
    set(h(i), 'Color', colors(i, :));  
    
    legendStrings2{i} = sprintf('$\\omega = %.2f$', omega_values(i)); 
    hold on;
end

legend(legendStrings2, 'Interpreter', 'latex', 'Location', 'best');

title('Error of $\phi$ vs $\phi_{linearized}$, being: $(\frac{\phi - \delta + \pi / 2}{\phi})100\% $ for t = [0, 3T]');
xlabel('time (t)');

ylabel('\% Error of $\phi$ vs $\phi_{linearized}$');

grid on;

% print(7, '-dpng', '-r600');

fig72 = figure(72);
set(fig72, 'Position', [360, 420, 700, 500]);
for i = 1:3
    tspan = linspace(0, 20*2* pi /omega_values(i),1000);

    Opt = odeset('Stats','on', 'MaxStep', 1e-2, 'Refine', 10, 'RelTol', 1e-12,'AbsTol', 1e-14);
    
    hold on;
    [t_x, X] = ode45(@(t, X) diffeq(t, X, g, omega_values(i), A, R), tspan, X0, Opt);
    phi  = X(:,1);
  
    [t, Y] = ode45(@(t, Y) diffeq1(t, Y, g, omega_values(i), A, R), tspan, Y0, Opt);
    delta = Y(:,1);

    relative_diff = (phi - delta + (pi/2)) .* 100 ./ phi;
    

    h(i) = plot(t, relative_diff, 'k'); 
    set(h(i), 'Color', colors(i, :));  
    
    legendStrings2{i} = sprintf('$\\omega = %.2f$', omega_values(i)); 
    hold on;
end

legend(legendStrings2, 'Interpreter', 'latex', 'Location', 'best');

title('Error of $\phi$ vs $\phi_{linearized}$, being: $(\frac{\phi - \delta + \pi / 2}{\phi})100\% $ for t = [0, 20T]');
xlabel('time (t)');

ylabel('\% Error of $\phi$ vs $\phi_{linearized}$');

grid on;

xlim([0, 10]);

% print(72, '-dpng', '-r600');



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



function Xdot = diffeq(t, X, g, omega, A, R)
    phi = X(1);
    phidot = X(2);
    Xdot = [phidot; ((-g*cos(phi))/R) - (A^2) * (omega^2) * (1 +cos(phi)) * (sin(omega*t)^2)*sin(phi)];
   

end

function Ydot = diffeq1(t, Y, g, omega, A, R)
    delta = Y(1);
    deltadot = Y(2);
    Ydot = [deltadot; (A^2 * omega^2 *(1-cos(2*omega*t))) / 2 - (g*delta)/R];
end