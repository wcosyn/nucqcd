# A(e,e'NN) with correlated pairs: unfinished since 2012, never used for results (docs/HISTORY.md).
# The classes were moved here from the Knockout (DoubleNModel, DoubleNCross) and MePhys (Pair,
# Newcoef) libraries in 2026-11; RecMosh (Moshinsky brackets) stays in MePhys.
nucqcd_add_program(twonucl SOURCES twonucl.cpp DoubleNModel.cpp DoubleNCross.cpp pair.cpp newcoef.cpp
                   LINK Knockout)
