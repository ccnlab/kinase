# Consolidated TODO list for Kinase

- [ ] This is the `LrnM` param in neuron-level `SpkCa` params.  However, using `MTau=10` `PTau=40` instead of the reverse could obviate the need for this.  TODO: test this more systematically, get rid of LrnM if possible, though used for optimization thresholding -- find a better way to do that?

- [ ] PThrMin is promising -- test more extensively.

- [ ] test SynSpkCa vs. SynContCa -- dwts are fairly correlated so likely similar but worth checking.  SynSpk diverges from NeurSpkCa pretty dramatically, though overall variance is similar.

- [ ] `SynNMDACa` should be using Urakbuo ITau=100, Tau=30 dynamics for *sender* NMDA used in learning!  This is not currently the case!  Snmda added, ITau is bad, but what about Tau?

- [ ] consider some combination of function of overall CaMKII level, and the error-driven diff -- e.g. multiply two factors?  something.. basically need to figure out how to get the STDP result plus the err diff, and perhaps have some additional mechanism of differentiation -- the biophysically grounded mech is pretty good overall, but CHL `NeurSpkCa` does better for avoiding hogging..

