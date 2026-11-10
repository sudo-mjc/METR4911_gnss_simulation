function m = eval_rmse(t, p_est, p_true, opts)
%EVAL_RMSE  Horizontal position accuracy metrics for the LC INS/GNSS sim.
%
%   m = eval_rmse(t, p_est, p_true)
%   m = eval_rmse(t, p_est, p_true, opts)
%
%   Checks the project goal:
%       ">= 85% of solutions within <= 0.2 m (RMSE) when pathing"
%
%   INPUTS
%     t       N x 1  time vector [s]
%     p_est   N x 2  estimated position [E N] in local frame [m]
%     p_true  N x 2  truth position     [E N] in same frame  [m]
%     opts    (optional struct)
%       .thr       error threshold [m]                       (default 0.2)
%       .goal      required pass fraction                    (default 0.85)
%       .win       rolling RMSE window length [s]            (default 1.0)
%       .mask      N x 1 logical, epochs to evaluate          (default all)
%                  e.g. gnss_valid -> GNSS-epoch accuracy,
%                       ~gnss_valid -> INS coast accuracy
%       .t_skip    initial convergence time to ignore [s]    (default 0)
%       .plot      true to plot error and rolling RMSE       (default true)
%
%   OUTPUT struct m
%     .rmse           overall horizontal RMSE [m]
%     .rmse_E/.rmse_N per-axis RMSE [m]
%     .err            N x 1 horizontal error per epoch [m] (NaN if excluded)
%     .pct_within     % of evaluated epochs with |error| <= thr
%     .rmse_win       N x 1 rolling-window RMSE [m]
%     .pct_win_within % of evaluated epochs with rolling RMSE <= thr
%     .cep50/.cep95   50th / 95th percentile horizontal error [m]
%     .max_err        worst horizontal error [m]
%     .pass           true if pct_win_within >= goal*100
%
%   Example (Simulink logged outputs):
%     t      = out.p_est.Time;
%     p_est  = squeeze(out.p_est.Data);          % make N x 2
%     p_true = interp1(truth.t, truth.p(:,1:2), t);  % align rates
%     m = eval_rmse(t, p_est, p_true, struct('t_skip', 10));

    if nargin < 4, opts = struct(); end
    thr    = getopt(opts, 'thr',    1.0);
    goal   = getopt(opts, 'goal',   0.85);
    win    = getopt(opts, 'win',    1.0);
    mask   = getopt(opts, 'mask',   true(numel(t),1));
    t_skip = getopt(opts, 't_skip', 0);
    doPlot = getopt(opts, 'plot',   true);

    % --- shape checks ---------------------------------------------------
    t = t(:);
    if size(p_est,1)  ~= numel(t), p_est  = p_est.';  end
    if size(p_true,1) ~= numel(t), p_true = p_true.'; end
    assert(size(p_est,1)==numel(t) && size(p_true,1)==numel(t), ...
        'p_est and p_true must have one row per time sample.');
    p_est = p_est(:,1:2);  p_true = p_true(:,1:2);

    % --- evaluation set: user mask, after convergence, finite values ----
    use = logical(mask(:)) & (t >= t(1) + t_skip) & ...
          all(isfinite(p_est),2) & all(isfinite(p_true),2);

    % --- per-epoch horizontal error -------------------------------------
    d   = p_est - p_true;              % [dE dN]
    e   = hypot(d(:,1), d(:,2));       % horizontal error magnitude
    eu  = e(use);

    m.rmse   = sqrt(mean(eu.^2));
    m.rmse_E = sqrt(mean(d(use,1).^2));
    m.rmse_N = sqrt(mean(d(use,2).^2));
    m.err    = nan(size(e));  m.err(use) = e(use);

    % --- metric 1: % of individual solutions within thr ------------------
    m.pct_within = 100 * mean(eu <= thr);

    % --- metric 2: rolling RMSE over a time window -----------------------
    % RMSE is a statistic over many samples, so "a solution within 0.2 m
    % RMSE" is best read as: the RMSE over a short window around it is
    % <= thr. Window is in seconds so it works for any sample rate.
    dt  = median(diff(t));
    nW  = max(1, round(win / dt));
    e2  = e.^2;  e2(~use) = NaN;
    m.rmse_win       = sqrt(movmean(e2, nW, 'omitnan'));
    m.pct_win_within = 100 * mean(m.rmse_win(use) <= thr);

    % --- distribution stats ---------------------------------------------
    m.cep50   = prctile(eu, 50);
    m.cep95   = prctile(eu, 95);
    m.max_err = max(eu);
    m.n_eval  = nnz(use);
    m.pass    = m.pct_win_within >= goal*100;

    % --- report ----------------------------------------------------------
    fprintf('\n--- Position accuracy (%d epochs, thr = %.2f m) ---\n', m.n_eval, thr);
    fprintf('Overall RMSE         : %.3f m  (E %.3f, N %.3f)\n', m.rmse, m.rmse_E, m.rmse_N);
    fprintf('CEP50 / CEP95 / max  : %.3f / %.3f / %.3f m\n', m.cep50, m.cep95, m.max_err);
    fprintf('Epochs |e| <= thr    : %.1f %%\n', m.pct_within);
    fprintf('Epochs %.1fs-RMSE<=thr: %.1f %%\n', win, m.pct_win_within);
    fprintf('Goal (>= %.0f %%)      : %s\n\n', goal*100, ternary(m.pass,'PASS','FAIL'));

    if doPlot
        figure('Name','Position accuracy');
        tiledlayout(2,1);
        nexttile; hold on; grid on;
        plot(t, m.err, 'DisplayName','|e| horizontal');
        plot(t, m.rmse_win, 'LineWidth',1.5, ...
             'DisplayName',sprintf('%.1f s rolling RMSE',win));
        yline(thr,'r--','DisplayName',sprintf('%.2f m threshold',thr));
        xlabel('Time [s]'); ylabel('Error [m]'); legend('Location','best');
        title(sprintf('RMSE %.3f m | %.1f%% within %.2f m', ...
              m.rmse, m.pct_win_within, thr));
        nexttile; grid on;
        histogram(eu, 50, 'Normalization','cdf');
        xline(thr,'r--'); yline(goal,'k:');
        xlabel('Horizontal error [m]'); ylabel('CDF');
        title('Error CDF (goal: curve above dotted line at threshold)');
    end
end

function v = getopt(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function o = ternary(c, a, b)
    if c, o = a; else, o = b; end
end