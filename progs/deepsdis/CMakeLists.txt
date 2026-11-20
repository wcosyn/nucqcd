nucqcd_add_programs(LINK DIS MePhys PROGRAMS
  deepstest               # DEEPS-like reduced structure functions
  deepsfit                # fit of the FSI parameters to the DEEPS data
  bonustest               # BONuS-like ratios (numerator and denominator separately), with the data
                          # renormalized by the BONuS MC; used for the plots in the BONuS paper.
                          # Use this one, not bonusextrapolate(2)
  bonusscatterfit         # refit of the rescattering parameters, no sensible results (spectator momenta too low)
  bonusfitreadout         # writes out some of the data renormalization fit values
  bonusfit_rig            # idem
  bonusextrapolate        # BONuS-like F2n ratios: data points for fitscript.py
  bonusextrapolate2       # BONuS-like F2n ratios: extrapolation procedure in fitscript.py
  bonus_F2nfit            # pole extrapolation of the BONuS data (renormalized), also independent fit of all angles
  bonus_F2nfit_sys        # idem, used for the systematic errors
  bonus_F2nfit_allangles  # idem, fitting all angles at the same time: NOT USED, no good results
  hermes_comp             # HERMES b1 kinematics and possible mixing of structure functions
  b1_calc                 # VNA convolution b1, comparison with Shunzo Kumano's model
  b1_calc_comp            # idem
  band_calc)
