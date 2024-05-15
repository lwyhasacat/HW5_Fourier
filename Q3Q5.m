a_x = 1.0;

N = 200; 
dx = 0.0050; 
dt = 0.0006; 
mus = [0.1, 0.01, 0.0001]; 
frames = 500; 
tol = 1e-6; 

U = zeros(N+2, 1);
S = ones(N+2, 1); 

figure;
Us = repmat({zeros(N+2, 1)}, 1, length(mus));
for idx = 1:length(mus)
    subplot(1, length(mus), idx);
    caxs(idx) = plot(Us{idx});
    ylim([-1, 2]); 
    title(sprintf('μ = %.4f', mus(idx)));
end

video_filename = 'advection_diffusion_simulation.mp4';
video = VideoWriter(video_filename, 'MPEG-4');
open(video);

for frame = 1:frames
    steady_state = true;
    for idx = 1:length(mus)
        U_prev = Us{idx};
        Us{idx} = forward_euler_1d_tri(Us{idx}, dx, dt, mus(idx), S);
        set(caxs(idx), 'YData', Us{idx});
        % if steady state
        if max(abs(Us{idx} - U_prev)) > tol
            steady_state = false;
        end
    end
    drawnow;
    frame_data = getframe(gcf);
    writeVideo(video, frame_data);
    if steady_state
        break;
    end
end

close(video);

function U_next = forward_euler_1d_tri(U, dx, dt, mu, S)
    nx = length(U) - 2;
    U_next = U;
    F = U(2:end-1) + dt * S(2:end-1);
    U_next(2:end-1) = tridiagonal_solver(mu, dt, dx, F, nx);
end

% tridiagonal solver
function V = tridiagonal_solver(mu, dt, dx, F, nx)
    c = mu * (dt / dx^2);
    ab = zeros(nx, 3);
    ab(2:end, 1) = -c;
    ab(:, 2) = 1 + 2*c;
    ab(1:end-1, 3) = -c;
    V = spdiags(ab, -1:1, nx, nx) \ F(:);
end
