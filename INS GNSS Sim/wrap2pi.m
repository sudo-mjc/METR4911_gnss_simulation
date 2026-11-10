function psi_wrapped = wrap2pi(psi)
%WRAP2PI Wraps angle to range of [-pi pi]
% 
% Created to call wherever needed and wrap heading values or an array of values
%

    psi_wrapped = atan2(sin(psi), cos(psi));

end