# CaDP
# Change Log
All notable changes to this project will be documented in this file.

## 1.4 - 
### TO-DO
- The framework for Figure 4-5 is already finished. Update pdf with these results.
- Why is 3C not fully replicated? Figure this out.


## 1.3 - 2017-03-24
### Changed
- CaDP.m: Figure 3C is added with partial success. The shape of the STDP curves are similar but the amplitudes are not.
- CaDP.m: Figure 5A-B are succesfully added.
- Minor changes in the other files

### Added

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
