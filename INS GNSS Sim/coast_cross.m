function [coast_end, crossed, t_cross, T_coast, sig_start, sig_peak, m_hat] = ...
                                   coast_cross(sigma_p, gnss_on, gnss_valid, cross_prm)
%COAST_CROSS  Measures each coast: did sigma_p reach sigma_set, when, and
%                how fast was it growing?
%   A coast runs from receiver switch-off to the first valid fix after waking.
%   Its clock starts at 0 on the switch-off tick.
%     crossed   : sigma_p reached sigma_set during the coast
%     t_cross   : first time sigma_p >= sigma_set   [s]
%     T_coast   : time of the first fix             [s]
%     sig_start : sigma_p at switch-off             [m]
%     sig_peak  : highest sigma_p in the coast      [m]
%     m_hat     : measured growth exponent of D = sigma_p^2 - sig_start^2
%                 near the end of the coast, clamped to [2, 6].
%                 0 = no estimate (coast < 2 s or D not yet positive);
%                 sleep_ctrl then falls back to ctrl_prm.m_fix.
%   Outputs update on the coast_end tick and are held until the next one.
%
%   How m_hat is measured
%     sigma_p comes from P, which has no measurement noise in it, so two
%     points on the curve give the exponent exactly for that coast:
%         m_hat = ln(D_end / D_A) / ln(T_end / t_A)
%     The coast length is not known in advance, so D is sampled at doubling
%     checkpoints (1, 2, 4, 8, ... s) and the last two are kept. At the
%     end of the coast the older one (t_A) lies in [T_end/4, T_end/2), so
%     ln(T_end/t_A) >= ln 2 (well-conditioned) and the slope describes the
%     LATE part of the coast, which is the part being extrapolated.
%     Clamp [2, 6]: 2 = velocity-dominated growth, 6 = heading error through
%     f_b (gyro bias). Values outside that come from a stop's accel burst or
%     a cross-covariance transient, not a growth law to extrapolate from.
%#codegen
persistent was_on active k tx pk n_done cr_h tx_h Tc_h pk_h s0 s0_h
persistent tA DA tB DB tck m_h
if isempty(was_on)
    was_on = true;  active = false;  k = 0;  tx = -1;  pk = 0;  n_done = 0;
    cr_h = false;   tx_h = 0;  Tc_h = 0;  pk_h = 0; s0 = 0; s0_h = 0;
    tA = 0;  DA = 0;  tB = 0;  DB = 0;  tck = 1;  m_h = 0;
end
coast_end = false;

if was_on && ~gnss_on                 % receiver just switched off: coast starts
    active = true;  k = 0;  tx = -1;  pk = sigma_p;
    s0 = sigma_p;
    tA = 0;  DA = 0;  tB = 0;  DB = 0;  tck = 1;     % first checkpoint at 1 s
elseif active
    if gnss_valid                     % first fix: coast over
        active = false;  n_done = n_done + 1;
        cr_h = (tx >= 0);
        tx_h = max(tx, 0);
        Tc_h = (k + 1)*cross_prm.dt;        % this tick is k+1 ticks after switch-off
        pk_h = pk;
        s0_h = s0;

        % Growth exponent from the older checkpoint to the end of the coast
        D_end = pk^2 - s0^2;
        T_end = k*cross_prm.dt;             % last coast tick (pk includes it)
        if DA > 0 && D_end > DA && tA > 0
            m_h = log(D_end/DA) / log(T_end/tA);
            m_h = min(max(m_h, 2), 6);
        else
            m_h = 0;
        end

        coast_end = (n_done > cross_prm.n_skip);
    else
        k  = k + 1;
        pk = max(pk, sigma_p);
        if tx < 0 && sigma_p >= cross_prm.sigma_set
            tx = k*cross_prm.dt;            % first crossing of the target
        end
        if k*cross_prm.dt >= tck            % doubling checkpoint reached
            tA = tB;  DA = DB;
            tB = k*cross_prm.dt;  DB = sigma_p^2 - s0^2;
            tck = 2*tck;
        end
    end
end
was_on = gnss_on;

crossed = cr_h;  t_cross = tx_h;  T_coast = Tc_h;
sig_peak = pk_h; sig_start = s0_h; m_hat = m_h;

end