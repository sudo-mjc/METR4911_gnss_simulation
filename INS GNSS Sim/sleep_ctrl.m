function [T_sleep_cmd, dT] = sleep_ctrl(coast_end, crossed, t_cross, T_coast, ...
                                           sig_start, sig_peak, m_hat, ctrl_prm)
%SLEEP_CTRL  Updates the sleep length once per coast.
%   Error e = (ideal coast length) - (actual coast length):
%     Crossed     : e = t_cross - T_coast   (<= 0: the coast ran too long)
%     Not crossed : e = t_hat   - T_coast   (>= 0: room for a longer sleep)
%                   t_hat = T_coast * ((sigma_set^2 - sig_start^2)/D_end)^(1/m)
%   PI (adapt_en = true):  T = T_sleep0 + Kp*e + Ki*sum(e), clamped
%   Bang-bang (adapt_en = false): T = T -/+ dT_default, clamped (unchanged)
%
%   Changes from v1 (PI branch only; bang-bang baseline is untouched):
%   1. m from the data: m = m_hat (measured by coast_cross_v2) when
%      available, ctrl_prm.m_fix only as the fallback (m_hat = 0).
%   2. Sign of the proportional term: +Kp*e (v1 had -Kp*e, which pushed
%      against the integral term on every new error).
%   3. Anti-windup (conditional integration): the error is not added to
%      the integrator if that would push T further into a clamp it is
%      already against. Without it, a long outage (one coast of hundreds
%      of seconds -> large negative e) or a run of degraded coasts winds I
%      far negative, and T sits at T_min long after the fixes recover.
%   4. Unachievable coasts: if sig_start >= sigma_set, the fix itself was
%      already worse than the target (heavy degradation), so no sleep
%      length can meet it. T is commanded to T_min (keep sigma_p as low as
%      possible) but I is NOT updated, so when fixes recover the
%      controller returns to the sleep length it had learned before.
%#codegen
persistent T I dT_h
if isempty(T)
    T = ctrl_prm.T_sleep0; I = 0;  dT_h = 0;
end

if ctrl_prm.open_loop
    T = ctrl_prm.T_fix;
elseif coast_end
    T_prev = T;
    if ctrl_prm.adapt_en
        % ======== PI CONTROLLER ========
        if sig_start >= ctrl_prm.sigma_set
            % (4) target unachievable on this coast: be cautious, don't learn
            T = ctrl_prm.T_min;
        else
            if crossed
                e = t_cross - T_coast;
            else
                % (1) growth exponent: measured if available, else fixed
                m = ctrl_prm.m_fix;
                if m_hat > 0
                    m = m_hat;
                end
                D_end = sig_peak^2 - sig_start^2;
                if D_end > 0
                    t_hat = T_coast * ((ctrl_prm.sigma_set^2 - sig_start^2) ...
                                       / D_end)^(1/m);
                    e = t_hat - T_coast;
                else
                    e = 0;
                end
            end

            % (3) conditional integration: skip if it deepens saturation
            T_try  = ctrl_prm.T_sleep0 + ctrl_prm.Kp*e + ctrl_prm.Ki*(I + e);
            windup = (T_try > ctrl_prm.T_max && e > 0) || ...
                     (T_try < ctrl_prm.T_min && e < 0);
            if ~windup
                I = I + e;
            end

            % (2) correct sign on the proportional term
            T_raw = ctrl_prm.T_sleep0 + ctrl_prm.Kp*e + ctrl_prm.Ki*I;
            T = min(max(T_raw, ctrl_prm.T_min), ctrl_prm.T_max);
        end
    else
        % ======== BANG-bang Baseline (unchanged) ========
        if crossed
            T = T - ctrl_prm.dT_default;
        else
            T = T + ctrl_prm.dT_default;
        end
        T = min(max(T, ctrl_prm.T_min), ctrl_prm.T_max);
    end
    dT_h = T - T_prev;
end

T_sleep_cmd = T;  dT = dT_h;
end