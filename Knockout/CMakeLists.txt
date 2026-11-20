# libKnockout: quasi-elastic A(e,e'p), 2N knockout, rho/pion transparency, neutrino QE.
nucqcd_add_library(Knockout
  SOURCES
    DeuteronMomDistr.cpp
    Model.cpp
    Cross.cpp
    DoubleNModel.cpp
    DoubleNCross.cpp
    RhoTCross.cpp
    PionTCross.cpp
    RhoDeuteron.cpp
    WeakQEHadronCurrent.cpp
    WeakQECross.cpp
    DQEinclusive.cpp
  LINK MePhys Integrator Glauber)
