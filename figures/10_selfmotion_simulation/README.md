# Self-motion demo in simulation

The two self-motion experiments (`orbit-selfmotion-hold`, `orbit-selfmotion-moving`) in simulation: the end
effector starts on the inspection orbit (xi0 = 0) with the arm off the lift by beta0 = 0.5 rad, and the controller
returns the arm to the lift.

The simulation uses the exact model (the joints do exactly the commanded acceleration, as in the preflight), with
no friction and the paper's PD law. Compared with the robot experiment:
- no integral action;
- no fade-in of the controller (the robot's 1 s ramp);
- transverse poles at -8, beta poles at -4.

Result, both variants:
- beta follows the design beta0 (1 + 4t) e^{-4t} to within 4e-4 rad;
- the joints are back on the lift by about 2 s;
- the end effector stays on the path to within 1.2 um and 0.4 arcsec (|xi| <= 1.2e-6).

Regenerate (from FallTests/):

    g++ -O2 -std=c++17 -Iinclude -Iapps/robot -I/usr/include/eigen3 figures/10_selfmotion_simulation/simdemo.cpp \
        build/libse3pf.a -o /tmp/simdemo
    /tmp/simdemo figures/10_selfmotion_simulation
    python3 tools/plot_self_motion.py figures/10_selfmotion_simulation/orbit-selfmotion-hold_sim.csv --design 4 --sim
    python3 tools/plot_self_motion.py figures/10_selfmotion_simulation/orbit-selfmotion-moving_sim.csv --design 4 --sim

`*_sim.csv.preflight.csv` are links to the logs themselves: the plot takes the lift from the samples with beta = 0.
