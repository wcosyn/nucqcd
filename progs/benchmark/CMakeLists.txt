# Benchmarks (not tests): integrators compares rombergerN, cube_romb, cube_adaptive and Cuhre on test
# and real integrands (roadmap step 6). Results are recorded in docs/integrator_benchmark.md.
nucqcd_add_program(integrators LINK Knockout Glauber MePhys Integrator cuba)
