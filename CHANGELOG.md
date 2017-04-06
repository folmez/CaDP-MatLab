# CaDP
# Change Log
All notable changes to this project will be documented in this file.

## 1.9 - 2017-04-05
### Changed
- CaDP.m: Stochastic CaDP option is added to replicate Figure 2 in Shouval 2004.
- Currently running on INS server: CaDP('sCDAP_VST_det_and_sto', [-100:10:-20 -19:4:45 50:10:100], 20, 'save_workspace', 1);

### Added
- find_sCDAP_G_NMDA_gamma_dist_variance.m: If I am understanding correctly, variance of the gamma distribution, from which G_NMDA is randomly drawn, is determined so that the resulting Ca peak coffecient of variation is the same as the theory (Figure 1D in Shouval 2004).  This routine determines the variance so that the resulting Ca peak coeficient of variation is accurate with a relative error smaller than 0.01.
- display_sim_progress.m: A little routine for displaying simulation progress

## 1.8 - 2017-03-31
### Changed
- CaDP.m: Stochastic simulations tested. Results are not good.
- "stochastic" parameter option is added to "NMDAr_calcium_current.m".

### Added
- "generate_random_maximal_chord_conductance.m": As described in Shouval 2004.

## 1.7 - 2017-03-30
### Changed
- "Stochastic" parameter options are added to all reuqired routines except for "NMDAr_calcium_current.m".
- "Deterministic with stochastic parameters" parameter option is added to "NMDAr_calcium_current.m" to replicate Figure 1B in Shouval 2004. "Stochastic" version will be added later.

## 1.6 - 2017-03-29
### Changed
- CaDP.m: Figure S9 can also be replicated with the option 'VST' (varying spike times).
- Omega_calcium.m: "Stochastic" parameter options from Shouval 2004 are added.
- eta_calcium.m: "Stochastic" parameter options from Shouval 2004 are added.

## 1.5 - 2017-03-28
### Changed
- CaDP.m: It turns out that choice of dt=1 is a bit too crude. Yesterday's problems with Figure 4A is solved by choosing dt=0.1. When the step-size is too big, the BPAP decays too fast. So fast that 
- CaDP.m: Figure 4A is currently running on the INS cluster.
- CaDP.m: Figure 4B is partially replicated. For 1Hz, my output seems to be different that what's in the paper. But their result is also different from the Figure 3C. So they must be using different parameters without telling. So, I think, I am fine here.
- CaDP.m: Save workscape option is added to some choices to run on the INS cluster.

## 1.4 - 2017-03-27
### Changed
- CaDP.m: Figure 5C is succesfully added.
- Minor changes in the other files

## 1.3 - 2017-03-24
### Changed
- CaDP.m: Figure 3C is added with partial success. The shape of the STDP curves are similar but the amplitudes are not.
- CaDP.m: Figure 5A-B are succesfully added.
- Minor changes in the other files

## 1.2 - 2017-03-23
### Changed
- CaDP.m: Figure 3B is corrected with success. The problem was that I was considering half of all channels opening up after each spike. However, according to Shouval 2002, half of all closed channels open up after each spike. This means, for higher frequencies, while a portion of the channels are open, half of the remaining closed channels open up after a spike and this may make fraction of open NMDAr channels more than 0.5. This is why previously smaller frequencies produced good results but higher frequencies did not.

- NMDAr_calcium_current.m: A second option is added to calculate the number of open NMDAr channels according to above description. For every pre spike time, the number of closed channels right before the spike is recorded and used to track the number of channels that stay open due to a particular spike.

## 1.1 - 2017-03-22
### Changed
- CaDP.m: Figure 3B is added with partial success. The shape is similar but it is not quantitatively identical to what is in Shouval 2002. I will leave this as it is right now and continue with Figure 3C next time.
- calculate_spike_times.m: If the time interval is not long enough to the stimulation, simulation end time is updated.


## 1.1 - 2017-03-21
### Changed
- CaDP.m: Figure 3A is succesfully added
- update_w.m: Eta function has the units 1/sec rather than 1/msec. The synaptic plasticity equation is corrected.
- BPAP.m: "Narrow BPAP" and "BPAP + ADP" options are added according to Shouval 2002.
- Minor changes in the other files

### Added
- calculate_spike_times.m: This calculates presynaptic spike times in deterministic way.

## 1.0 - 2017-03-20
### Added
- CaDP.m: Figure 2 replication is succesfully added.
- NMDAr_calcium_current.m
- BPAP.m
- EPSP.m
- update_w.m
- update_Ca.m
- eta_calcium.m
- Omega_calcium.m
