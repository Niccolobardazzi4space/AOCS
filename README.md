# AOCS

## Detumbling phase
During this phase the spacecraft, just injected into the orbit had to be firstly slowed down till zero angular velocity.

## Slew manoeuvre
The spacecraft had to compute a manoeuvre to be aligned with the prescribed attitude, which in this case makes the spacecraft pointing the Earth.

## Earth pointing
In spite of the perturbations added which are SRP, magnetic field, gravity gradient and drag, the spacecraft had to preserve its motion pointing the Earth with the best accuracy as possible. For this reason it has been used a non-linear filter, the Extended Kalman Filter which provides very high performances. However a determination process is required at the beginning before starting the non-linear estimation.

### Installation
The whole project has been programmed in MATLAB language using Simulink.
