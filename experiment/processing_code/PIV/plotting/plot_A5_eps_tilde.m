function plot_A5_eps_tilde(data)
% A.5  Characteristic steepness eps_tilde(t) = sqrt(2 <eta_x^2>(t)).
% Wagner et al. 2023 Eq. 2.1 reproduction.

eta   = data.eta;
x_eta = data.x_eta;
t_eta = data.t_eta;
Fs    = data.Fs;
dx    = mean(diff(x_eta));

[Nx, Nt] = size(eta);
eps_t = zeros(Nt, 1);
for it = 1:Nt
    eta_x = gradient(eta(:, it), dx);
    eps_t(it) = sqrt(2 * mean(eta_x.^2, 'omitnan'));
end

% Smooth ~1 s
nroll = max(1, round(Fs));
eps_s = movmean(eps_t, nroll, 'omitnan');

fig = figure('Name','A5 eps_tilde (PIV)','Position',[100 100 1000 500],'Color','w');
plot(t_eta, eps_s, 'b-', 'LineWidth', 1.6);
xlabel('t (s)'); ylabel('\epsilon\sim(t)');
title('A.5  Characteristic steepness  \epsilon\sim = \surd(2 <\eta_x^2>)');
grid on;
save_figure(fig, 'A5_eps_tilde');
end
