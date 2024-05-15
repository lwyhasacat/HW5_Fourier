% Credit: Received help from Christine
clc, clear, close all;

L = 1;
t_final = 0.5;
a = 1;
S = 1;
mu_values = [0.1, 0.01];
dx_values = [0.1, 0.05, 0.01];

% plotting numerical solutions
figure(1);
sgtitle('Numerical Solutions at t = 0.5');

for mu_idx = 1:length(mu_values)
    mu = mu_values(mu_idx);
    for dx_idx = 1:length(dx_values)
        dx = dx_values(dx_idx);
        subplot(length(mu_values), length(dx_values), (mu_idx-1)*length(dx_values) + dx_idx);
        plot_numerical_solution(L, dx, mu, S, t_final, a);
    end
end

% plotting m(theta)
figure(2);
sgtitle('|1+\Delta t*m(\theta)| for different \mu and \Delta x values');

plot_m_theta(a, mu_values, dx_values);

function U_next = forward_euler(U, dx, dt, mu, S, a)
    U_next = U;
    for j = 2:length(U)-1
        U_next(j) = U(j) + dt * (-a / (2 * dx) * (U(j+1) - U(j-1)) + ...
                                 mu / (dx^2) * (U(j+1) - 2*U(j) + U(j-1)) + S);
    end
end

function plot_numerical_solution(L, dx, mu, S, t_final, a)
    n = floor(L / dx);
    x = linspace(0, L, n+2);
    x_internal = x(2:end-1); 
    dt_diffusion = 0.5 * dx^2 / mu;
    dt_advection = dx / a;
    dt = min(dt_diffusion, dt_advection);
    num_steps = floor(t_final / dt);
    
    fprintf('mu = %.3f, dx = %.3f, dt = %.5f, steps = %d\n', mu, dx, dt, num_steps);
    fprintf('CFL: mu*dt/dx^2 = %.3f\n', mu * dt / dx^2);
    
    U = zeros(1, n + 2);
    
    for k = 1:num_steps
        U = forward_euler(U, dx, dt, mu, S, a);
    end
    
    plot(x, U, '-o');
    title(sprintf('\\mu = %.3f, \\Delta x = %.3f, \\Delta t = %.4f', mu, dx, dt));
    xlabel('x');
    ylabel('U(x,t)');
end

function plot_m_theta(a, mu_values, dx_values)
    theta = linspace(0, 2 * pi, 400);
    
    for i = 1:length(mu_values)
        for j = 1:length(dx_values)
            mu = mu_values(i);
            dx = dx_values(j);
            dt_diffusion = 0.5 * dx^2 / mu;
            dt_advection = dx / a;
            dt = min(dt_diffusion, dt_advection);
            % plotting |abs(m(theta))|
%             m_theta = abs(-1i * dt * (a / dx) * sin(theta) + dt * (mu / dx^2) * (2 * cos(theta) - 2));
%             max_m = max(m_theta);
%            
%             subplot(length(mu_values), length(dx_values), (i-1)*length(dx_values) + j);
%             plot(theta, m_theta, 'LineWidth', 2, 'DisplayName', sprintf('|m(\\theta)|, max=%.2f', max_m));
%             hold on;
%             yline(1, 'r--', 'DisplayName', 'Stability limit (|m(\theta)|=1)');
%             hold off;
%             title(sprintf('\\mu = %.3f, \\Delta x = %.3f, \\Delta t = %.5f', mu, dx, dt));
%             xlabel('\theta');
%             ylabel('|m(\theta)|');
%             legend('Location', 'southoutside');
            
            % plotting |abs(m(theta))|
            m_theta = -1i / dx * sin(theta) + mu / dx^2 * (2 * cos(theta) - 2);
            fac = abs(1+dt*max(m_theta));
            
            subplot(length(mu_values), length(dx_values), (i-1)*length(dx_values) + j);
            plot(theta, m_theta, 'LineWidth', 2, 'DisplayName', sprintf('|1+\\Delta t*m(\\theta)|, max=%.2f', fac));
            hold on;
            yline(1, 'r--', 'DisplayName', 'Stability limit (|1+\Delta t*m(\theta)|=1)');
            hold off;
            title(sprintf('\\mu = %.3f, \\Delta x = %.3f, \\Delta t = %.5f', mu, dx, dt));
            xlabel('\theta');
            ylabel('|1+\Delta t*m(\theta)|');
            legend('Location', 'southoutside');
        end
    end
    
    set(gcf, 'Position', [100, 100, 1000, 800]);
end
