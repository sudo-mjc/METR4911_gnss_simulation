# METR4911 – Adaptive GNSS Duty-Cycling for INSGNSS Pedestrian Navigation

Simulation files for my METR4911 thesis (University of Queensland).
A loosely coupled INSGNSS error-state Kalman filter is used to adaptively control
GNSS duty cycle based on positional uncertainty, reducing power use while
maintaining positioning accuracy.

## Requirements
- MATLAB  Simulink R2024b (Minimum)
- Sensor Fusion and Tracking Toolbox
- Navigation Toolbox
- Control System Toolbox

## Running the simulation
Open this folder in MATLAB and run `sim_main.m`. This runs the full simulation
end to end trajectory generation, sensor simulation, filtering, duty-cycle
control and results.

## Customisation
Simulation settings are in `sim_config.m`, including trajectory, sensor noise,
GNSS degradation (`cfg.gnss.deg.en`) and controller mode (`adapt_en` for
adaptive PI vs fixed bang-bang). Edit these, then rerun `sim_main.m`.

## Related
PCB design files https://github.com/sudo-mjc/Linear-Charger-Power-PCB

## Author
Mitchell Crawford – Bachelor of Engineering (Mechatronics), UQ
