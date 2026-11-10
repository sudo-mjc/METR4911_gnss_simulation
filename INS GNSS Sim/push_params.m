function push_params(cfg)
%PUSH_PARAMS  Rebuild the flat block parameters and put them in the base workspace.
    prm = build_sim_params(cfg);
    f = fieldnames(prm);
    for k = 1:numel(f), assignin('base', f{k}, prm.(f{k})); end
end