# Rho production on the deuteron (coherent and incoherent, t- and z-integrated): an abandoned calculation, never used for
# results. RhoDeuteron was moved here from the Knockout library, and rhodint from progs/RhoTransparency, in 2026-11.
# The nested Romberg routine they use is progs/attic/rombergerN.hpp (it left Tools/Core with this move).
nucqcd_add_program(rhodint SOURCES rhodint.cpp RhoDeuteron.cpp LINK Knockout Glauber MePhys)
