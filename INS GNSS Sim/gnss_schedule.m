function [gnss_on, gnss_valid, T_sleep_active] = gnss_schedule(T_sleep_cmd, wake_early, sched_prm)
%GNSS_SCHEDULE  Receiver state machine: WARMUP -> ACTIVE -> SLEEP -> REACQ -> ACTIVE ...
%   gnss_on    : receiver powered (energy / duty cycle)
%   gnss_valid : one-tick pulse per fix at f_gnss (drives the ESKF update)
%   T_sleep_cmd is latched only on entry to SLEEP.
%#codegen
persistent mode k n_fix T_lat
WARMUP = uint8(0); ACTIVE = uint8(1); SLEEP = uint8(2); REACQ = uint8(3);
if isempty(mode)
    mode = WARMUP; k = 0; n_fix = 0; T_lat = sched_prm.T_sleep0;
end

k_fix   = round(1/(sched_prm.f_gnss*sched_prm.dt));   % ticks between fixes
N_fix   = round(sched_prm.T_active*sched_prm.f_gnss); % fixes per active window
k_warm  = round(sched_prm.T_warmup/sched_prm.dt);
k_reacq = round(sched_prm.t_reacq/sched_prm.dt);

gnss_valid = false;
switch mode
    case WARMUP
        gnss_on    = true;
        gnss_valid = (mod(k, k_fix) == 0);
        k = k + 1;
        if k >= k_warm
            mode = ACTIVE; k = 0; n_fix = 0;
        end
    case ACTIVE
        gnss_on    = true;
        gnss_valid = (mod(k, k_fix) == 0);  % first fix on the first ACTIVE tick
        k = k + 1;
        if gnss_valid, n_fix = n_fix + 1; end
        if n_fix >= N_fix
            mode  = SLEEP; k = 0;
            T_lat = T_sleep_cmd;            % latch at the falling edge
        end
    case SLEEP
        gnss_on = false;
        k = k + 1;
        if k >= round(T_lat/sched_prm.dt) || wake_early
            mode = REACQ; k = 0;
        end
    otherwise   % REACQ: powered, acquiring, no fix yet
        gnss_on = true;
        k = k + 1;
        if k >= k_reacq
            mode = ACTIVE; k = 0; n_fix = 0;
        end
end
T_sleep_active = T_lat;
end