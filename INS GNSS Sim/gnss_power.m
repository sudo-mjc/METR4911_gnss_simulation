function pw = gnss_power(gnss_on, opts)
%GNSS_POWER  Receiver energy and average power from a logged gnss_on signal.
%   pw = gnss_power(out.gnss_on)
%   pw = gnss_power(out.gnss_on, t_reacq=1, V=3.3, verbose=false)
%
%   Model (u-blox ZED-F9P datasheet, Table 13 / section 4.2, GPS+GLO+GAL+BDS):
%     REACQ  : first t_reacq s after each wake      I_acq = 90 mA
%     ACTIVE : rest of powered time (tracking)      I_trk = 85 mA
%     SLEEP  : receiver off, software backup        I_slp = 1.4 mA
%   Each rising edge of gnss_on is one reacquisition. The initial warm-up
%   (gnss_on already high at t = 0) is charged at tracking power.
%
%   Inputs
%     gnss_on : timeseries (logical/0-1), e.g. out.gnss_on
%     opts    : name-value overrides, see the arguments block below
%   Output struct pw
%     E [J], P_avg [W], I_avg [A], Q [mAh], duty [-], nWake, T_run [s],
%     T_on [s], T_off [s], dt [s]

arguments
    gnss_on
    opts.V       (1,1) double = 3.3      % [V]  supply voltage
    opts.I_acq   (1,1) double = 130e-3    % [A]  acquisition current
    opts.I_trk   (1,1) double = 85e-3    % [A]  tracking current
    opts.I_slp   (1,1) double = 1.4e-3   % [A]  software-backup current
    opts.t_reacq (1,1) double = 2        % [s]  hot start (datasheet); match the sim's prm.t_reacq
    opts.verbose (1,1) logical = true    % print a one-line summary
end

on = logical(gnss_on.Data(:));
dt = mean(diff(gnss_on.Time));           % sample period of the log [s]

P_acq = opts.V*opts.I_acq;
P_trk = opts.V*opts.I_trk;
P_slp = opts.V*opts.I_slp;

nWake = sum(diff(on) == 1);              % rising edges = reacquisitions
T_on  = sum(on)*dt;
T_off = sum(~on)*dt;
T_run = numel(on)*dt;

% Powered time at tracking power, plus the extra cost of each reacquisition window
E = P_trk*T_on + (P_acq - P_trk)*nWake*opts.t_reacq + P_slp*T_off;   % [J]

pw.E     = E;
pw.P_avg = E/T_run;                      % [W]
pw.I_avg = pw.P_avg/opts.V;              % [A]
pw.Q     = pw.I_avg*T_run/3.6;           % [mAh]  (A*s / 3.6 = mAh)
pw.duty  = T_on/T_run;
pw.nWake = nWake;
pw.T_run = T_run;  pw.T_on = T_on;  pw.T_off = T_off;  pw.dt = dt;
    
    if opts.verbose
        fprintf(['E = %.1f J   P_avg = %.1f mW   I_avg = %.1f mA   Q = %.2f mAh   ' ...
                 'duty = %.1f %%   wakes = %d\n'], ...
                E, 1e3*pw.P_avg, 1e3*pw.I_avg, pw.Q, 100*pw.duty, nWake);
    end
end