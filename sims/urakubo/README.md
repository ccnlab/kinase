# Urakubo = Urakubo et al, 2008

This model re-implements the highly biophysically realistic model of spike-timing dependent plasticity from:

* Urakubo, H., Honda, M., Froemke, R. C., & Kuroda, S. (2008). Requirement of an allosteric kinetics of NMDA receptors for spike timing-dependent plasticity. The Journal of Neuroscience, 28(13), 3310–3323. http://www.ncbi.nlm.nih.gov/pubmed/18367598 | [Main PDF](https://github.com/emer/axon/blob/master/examples/urakubo/papers/UrakuboHondaFroemkeEtAl08.pdf) | [Supplemental PDF](https://github.com/emer/axon/blob/master/examples/urakubo/papers/UrakuboHondaFroemkeEtAl08_suppl.pdf) | [Model Details](https://github.com/emer/axon/blob/master/examples/urakubo/papers/UrakuboHondaFroemkeEtAl08_model_sup.pdf) (last one has all the key details but a few typos and the code is needed to get all the details right).

This model captures the complex postsynaptic chemical cascades that result in the changing levels of AMPA receptor numbers in the postsynaptic density (PSD) of a spine, as a function of the changing levels of calcium (Ca2+) entering the spine through NMDA and VGCC (voltage-gated calcium) channels.  The model mechanisms are strongly constrained bottom-up by the known kinetics of these chemical processes, and amazingly it captures many experimental results in the LTP / LTD experimental literature.  Thus, it provides an extremely valuable resource for understanding how plasticity works and how realistic but much simpler models can be constructed, which are actually usable in larger-scale simulations.

The complexity of this model of a single synaptic spine involves 116 different variables updated according to interacting differential equations requiring a time constant of 5e-5 (20 steps per the typical msec step used in axon).  In previous work with the original Genesis implementation of this model (see [CCN textbook](https://CompCogNeuro.org)), we were able to capture the net effect of this model on synaptic changes, at a rate-code level using random poisson spike trains, using the simple "checkmark" linearized version of the BCM learning rule (called "XCAL" for historical reasons), with an `r` value of around .9 (81% of the variance) -- see figure below.  However, with the actual spiking present in axon, there is an important remaining question as to whether a more realistic (yet still simplified) learning mechanism might be able to take advantage of some aspects of the detailed spiking patterns to do "better" somehow than the simplified rate code model.

The original Genesis model is a good object lesson in how a generic modeling framework can obscure what is actually being computed.  The model shows up as a complex soup of crossing lines in the Genesis GUI, and the code is a massive dump of inscrutable numbers, with key bits of actual computation carefully hidden across many different bits of code, with seemingly endless flags and conditions making it very difficult to determine what code path is actually in effect.  In addition, Genesis is ancient software that requires an old gcc compiler to build it, and it only runs under X windows.

By contrast, the present model is purpose-built in compact, well-documented highly readable Go code, using a few very simple types in the [emergent chem](https://github.com/emer/emergent/tree/master/chem) package, providing the basic `React` and `Enz` mechanisms used for simulating the chemical reactions.

# Model Summary

![Urakubo et al., (2008) Model](figs/fig_urakubo_et_al_08_bcm_sum.png?raw=true "Summary of Urakubo et al., 2008 model, and how it simulates many STDP data points, and behavior for more realistic spike trains can be summarized by a simple linearized BCM-like function")

The above figure summarizes the model and some of the STDP data that it simulates, along with the derivation of the XCAL checkmark linearized BCM-like learning function from more realistic spiking inputs.  The detailed CaMKII portion of the model is directly from Dupont et al (2003) with some parameter modifications, while the compartmental spiking model is from Migliore et al. (1999) and Poirazi et al. (2003).  Other signaling cascades were adapted from Bhalla and Iyengar (1999) and earlier cerebellar LTD models by this group: Kuroda et al. (2001); Doi et al. (2005).

In our version, we omitted the full compartmental model in favor of the basic Axon AdEx-style model which just summarizes the Vm curve resulting from a spike.  However, we did have to adjust a few parameters to match the Vm curve from the Genesis model (though these do not appear to make much difference in the end).

We also added VGCC channels (not otherwise used in axon -- verified that they do make a difference for the basic STDP results), and used the allosteric NMDA dynamics from the original Urakubo model instead of the simpler NMDA channel used in axon.  See the `NMDAAxon` option.

# Basic Usage

Toolbar at the top does actions (mouse over and hold to see tooltips), e.g.:
* `Init` re-initializes with any params updates
* `Run` runs current `Stim` program

Params list at left side control lots of things (see tooltips), e.g.:
* `Stim` determines the stimulation program -- see `stims.go` for code -- figures below show some examples.
* `ISISec` = spacing between individual stimulation patterns
* `NReps` = number of repetitions of patterns -- more = bigger effects on weights
* `FinalSecs` = number of seconds to allow AMPAR trafficking to settle out -- has significant effects on final results, particularly on showing LTD effects which emerge over time.  50 secs is long enough to be close to asymptotic -- 20 vs. 50 still shows some diffs but 50 vs. 100 is pretty close.
* Other params used in different cases -- see examples below -- and below that are overall parameters of the neuron etc.

Tab bar selects different plots:
* `DWtPlot` shows final delta-weights (change in synaptic weight), i.e., `Trp_AMPAR` = trapped AMPA receptors, as a proportion of `InitWt` which is baseline value -- only updated for the `Sweep` Stims.
* `PhaseDWtPlot` shows DWt for `ThetaErr` case, which manipulates minus and plus phase activity states according to the Contrastive Hebbian Learning (CHL) algorithm.
* `Msec*Plot` shows all the model state values evolving over time at different temporal scales (100 = 1 point for every 100 msec, etc)

See https://github.com/emer/axon/blob/master/examples/urakubo/results for plots and tab-separated-value (TSV) data for various cases, some of which are summarized below.

# Classical STDP Replication

The classic Bi & Poo STDP curve emerges reliably with e.g., 100 reps with 50 secs final "burn in" of the weights -- maximum burn-in occurs after about 100 secs and after then things can start to decay. The curve gets stronger with more reps as well.  The X axis is Tpost - Tpre, so the "causal" pre-post ordering is positive numbers on the right.

![STDPSweep](results/fig_urakubo22_stdp_100rep_final50.png?raw=true "Classic STDP with standard params -- X axis is Tpost - Tpre")

* `Stim = STDPSweep`
* `NReps = 100`
* `FinalSecs = 50`

# STDP in Packets

What happens in a slightly more "realistic" scenario where pre-post spikes occur at more naturalistic overall frequencies? This was tested with the STDPPacketSweep inputs: 200 msec of 25 or 50 hz spiking, with fixed pre-post time offsets.  This case represents a kind of "best case" for whether there might be some kind of STDP-like signal left in more realistic spike trains -- it is still highly unlikely that in the noisy environment of the brain, you would ever see such reliable pre-post timing across a sequence of spikes like this.

Interestingly, this shows the same qualitative effects but in a much weaker form -- this means that there may be *some* residual STDP effect in the overall learning rule -- can't quite conclude that yet because this model does not capture CHL.

![STDPPacketSweep](results/fig_urakubo22_stdpacket_25hz_isi1_10rep_final50.png?raw=true "STDP in theta-window packets of pre-post activity")

* `Stim = STDPPacketSweep`
* `NReps = 10`
* `FinalSecs = 50`
* `SendHz = 25 or 50`
* `DurMsec = 200`

# XCAL / BCM checkmark learning rule

The data used to derive the XCAL check-mark linearized version of the BCM curve manipulated send and recv Hz and also overall duration -- the main plots are for a given sending Hz, varying recv Hz (X axis) and duration.  These did NOT use poisson random spiking but rather ensured that pre and post fired 180 degrees out phase, so as to effectively negate any timing issues (equally pre or post).  Here are some plots for different levels of sending activity:

![XCAL Original](figs/fig_urakubo_xcal_3d_fits.png?raw=true "Original XCAL data from Genesis model (black underlying grid), with fitted checkmark function overlaid in red on top")

Here is the replication in current model -- viewed as stacked slices through the Z depth axis of the central panel of the above 3D plots -- the X axis is Recv Hz, different lines as shown in legend are  different durations of overall firing.

Sender = 100Hz:

![OppPhaseDur 100 SendHz](results/fig_urakubo22_opphasedur_100shz_final50_1rep.png?raw=true "Opposite phase firing (180 degrees out of phase at start) as function of Recv Hz on X axis and different durations (lines)")

Sender = 50Hz:

![OppPhaseDur 50 SendHz](results/fig_urakubo22_opphasedur_50shz_final50_1rep.png?raw=true "Opposite phase firing (180 degrees out of phase at start) as function of Recv Hz on X axis and different durations (lines)")

* `Stim = OpPhaseDurSweep`
* `NReps = 1` 
* `FinalSecs = 50`
* `SendHz = 50 or 100`

# Random Poisson Firing

This case shows what bad memory had consolidated as the basis for the original XCAL checkmark learning rule -- random poisson firing of send and recv neurons over a packet of a given duration. This produces a more complex curve, e.g., for 200msec activity packets (theta cycle), as a function of recv firing rate (X axis) and sending firing rate (legend / different lines):

![Poisson 200msec](results/fig_urakubo22_poisson_200msec_isi1_final50_10rep.png?raw=true "Poisson random 200 msec packets of pre - post activity at different rates -- X is Recv Hz, Send Hz are different lines as shown in legend")

* `Stim = PoissonHzSweep`
* `NReps = 10` 
* `FinalSecs = 50`
* `DurMsec = 200`

# CHL Error-driven Learning

The `ThetaErr` Stim runs the current Zito lab empirical experiment on the model, testing for plus-phase vs. minus phase differences in firing rate over a 200msec window (theta frequency), with the minus phase the first 100 msec, plus phase the second 200 msec.

![CHL Minus-Plus](figs/fig_chl_minus_plus_dwt_on_off.png?raw=true "Qualitative picture of the CHL error-driven learning function over time, as a minus-phase prediction followed by a plus-phase outcome.  On the Off-Off and On-On diagonal, the weights do not change, and on the off-diagonal terms you get either positive weight change (LTP, for a miss as Off->On) or negative weight change (LTD for an error as On -> Off)")

The X axis on the plot is CHL error signal for each condition: 

* CHL = X+Y+ - X-Y- 

With X = normalized frequency for sender (Hz/100), Y = norm recv frequency.  First term is hebbian co-product in the plus phase (correct answer, outcome, target), second is for the minus phase (guess, prediction, actual model output).

The Y axis is the DWt result of the model. If the model was perfectly capturing the CHL learning rule, the points would all be along a positive diagonal line (i.e., the r = 1 line).  Clearly, the original Urakubo model does not accurately capture the CHL function, although some subset of points are not too far off.  It is particularly inaccurate on the negative error cases, which require strong LTD, despite relatively high levels of activity.

![ThetaErr RSPk 200ms](results/old/fig_urakubo22_thetaerr_rspk_200ms_final100_10rep.png?raw=true "Theta (200 msec) packets of pre - post activity at different rates -- X is Recv Hz, Send Hz are different lines as shown in legend -- recv unit is driven by punctate spikes")

* `Stim = ThetaErr`
* `ISISec = 1`
* `NReps = 10` 
* `FinalSecs = 50`
* `DurMsec = 200`
* `RGClamp = false or true`

The next plot shows the case for `RGClamp = true`, where the recv unit is current-clamped to drive its spiking, instead of phasically driving individual spikes at a prescribed regular interval according to the target Hz.

![ThetaErr RGe 200ms](results/fig_urakubo22_thetaerr_rge_200ms_final50_10rep.png?raw=true "Theta (200 msec) packets of pre - post activity at different rates -- X is Recv Hz, Send Hz are different lines as shown in legend -- recv unit is driven by current clamp")

# Basic model behavior

You can see single stimulus events in the `MsecPlot` at full temporal resolution by running the `STDP` or `Poisson` Stim cases with `Nreps = 1, ISISec = 0, FinalSecs = 0`, so you just see the one trace.  This can be used for exploring parameter differences, seeing the Ca-driven dynamics etc.

# CaMKII drives everything?

![STDP Chem](results/fig_urakubo22_stdp_100rep_final1_chem.png?raw=true "Signaling state at end of std STDP induction")

The above figure suggests that the final dWt is strongly related to CaMKII levels. 

TODO: measure peak Ca too and plot that!  have to grab it within each spike window.


# LFS 900 @ 1Hz = LTD?

One standard way of generating LTD is low-frequency stimulation (LFS) for a large number of repetitions, e.g., 900 @ 1Hz (Dudek & Bear, 1992).  However, this does not produce the appropriate results in the basic Urakubo model:

![LFS 900 1hz](results/fig_urakubo22_lfs900_1hz_final20_ge.1.png?raw=true "LFS 900 reps at 1Hz, should cause LTD but does not")

* `Stim = STDP`
* `NReps = 900` 
* `FinalSecs = 1`
* `ISISec = 1`
* `DeltaT = 0` // key!
* `GeStim = .1 or .2` // .1 produces 2-5 mV EPSP, .2 = 10 mV

# References

* Akyol, Z., Bartos, J. A., Merrill, M. A., Faga, L. A., Jaren, O. R., Shea, M. A., & Hell, J. W. (2004). Apo-Calmodulin Binds with its C-terminal Domain to the N-Methyl-d-aspartate Receptor NR1 C0 Region *. Journal of Biological Chemistry, 279(3), 2166–2175. https://doi.org/10.1074/jbc.M302542200

* Barcomb, K., Hell, J. W., Benke, T. A., & Bayer, K. U. (2016). The CaMKII/GluN2B Protein Interaction Maintains Synaptic Strength. *Journal of Biological Chemistry, 291(31),* 16082–16089. https://doi.org/10.1074/jbc.M116.734822

* Bayer, K.-U., De Koninck, P., Leonard, A. S., Hell, J. W., & Schulman, H. (2001). Interaction with the NMDA receptor locks CaMKII in an active conformation. *Nature, 411(6839),* 801–805. https://doi.org/10.1038/35081080

* Bhalla, U. S., & Iyengar, R. (1999). Emergent properties of networks of biological signaling pathways. *Science.* https://doi.org/10.1126/science.283.5400.381

* Coultrap, S. J., & Bayer, K. U. (2012). CaMKII regulation in information processing and storage. *Trends in Neurosciences, 35(10),* 607–618. https://doi.org/10.1016/j.tins.2012.05.003

* Dudek, S. M., & Bear, M. F. (1992). Homosynaptic long-term depression in area CA1 of hippocampus and effects of N-methyl-D-aspartate receptor blockade. *Proceedings of the National Academy of Sciences U. S. A., 89(10),* 4363–4367. http://www.ncbi.nlm.nih.gov/pubmed/1350090

* Doi, T., Kuroda, S., Michikawa, T., & Kawato, M. (2005). Inositol 1,4,5-trisphosphate-dependent Ca2+ threshold dynamics detect spike timing in cerebellar purkinje cells. *Journal of Neuroscience, 25(4),* 950–961. https://doi.org/10.1523/JNEUROSCI.2727-04.2005

* Dupont, G., Houart, G., & De Koninck, P. (2003). Sensitivity of CaM kinase II to the frequency of Ca2+ oscillations: A simple model. *Cell Calcium, 34(6),* 485–497. https://doi.org/10.1016/S0143-4160(03)00152-0

* Kuroda, S., Schweighofer, N., & Kawato, M. (2001). Exploration of signal transduction pathways in cerebellar long-term depression by kinetic simulation. *Journal of Neuroscience, 21(15),* 5693–5702. https://doi.org/10.1523/JNEUROSCI.21-15-05693.2001

* Migliore, M., Hoffman, D. A., Magee, J. C., & Johnston, D. (1999). Role of an A-Type K+ Conductance in the Back-Propagation of Action Potentials in the Dendrites of Hippocampal Pyramidal Neurons. *Journal of Computational Neuroscience, 7(1),* 5–15. https://doi.org/10.1023/A:1008906225285

* Poirazi, P., Brannon, T., & Mel, B. W. (2003). Arithmetic of Subthreshold Synaptic Summation in a Model CA1 Pyramidal Cell. *Neuron, 37(6),* 977–987. https://doi.org/10.1016/S0896-6273(03)00148-X

* Urakubo, H., Honda, M., Froemke, R. C., & Kuroda, S. (2008). Requirement of an allosteric kinetics of NMDA receptors for spike timing-dependent plasticity. *The Journal of Neuroscience, 28(13),* 3310–3323. http://www.ncbi.nlm.nih.gov/pubmed/18367598

* Umemiya, M., Chen, N., Raymond, L. A., & Murphy, T. H. (2001). A Calcium-Dependent Feedback Mechanism Participates in Shaping Single NMDA Miniature EPSCs. Journal of Neuroscience, 21(1), 1–9. https://doi.org/10.1523/JNEUROSCI.21-01-00001.2001

