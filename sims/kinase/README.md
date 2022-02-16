# Kinase

This model is developing a new learning mechanisms based on the interactions between two kinases: CaMKII and DAPK1, which promote LTP and LTD respectively.   Computationally, DAPK1 represents the minus phase, and CaMKII represents the plus phase, with the subtractive relationship as expressed in the CHL learning algorithm:

```
    dW = (X^+ Y^+) - (X^- Y^-)
```

Where X = sending activation, Y = receiving, and superscripts indicate plus phase (outcome, target state) and minus phase (prediction, guess state) respectively.

The recent data from the Zito lab confirms the basic predictions of this learning rule, but many questions remain about the detailed underlying mechanisms and computational implications.

The model is being developed at multiple levels, with the two overarching goals of understanding biological mechanisms and exploring computational implications / advantages thereof.

* **KinaseBMax:** The most biophysically-detailed model is based on Urakubo et al, 2008, augmented with GluN2B binding of the two kinases, producing a competition between LTP and LTD.  The GluN2B binding for CaMKII has worked well, and DAPK1 is largely guesswork but currently is showing integration behavior, and a slower time constant than CaMKII, but it is not sufficiently bounded and it is unclear how to scale its dynamics relative to CaMKII.  Also, it is unclear how the other PP factors (PP1, PP2A, CaN) play into the story -- these might be important for scoping the learning window and altering dynamics.

* **KinaseBx:** Various potential levels of abstraction building up from from the detailed model, bridging to a high-level abstract model that can be run on large scale models.

* **Kinase:** The fully abstracted model that runs at scale, which will ultimately be what we call the "Kinase learning rule".  This will be a work in progress, starting with various levels of abstraction, the first of which is **KinaseAMax** that is the most abstract / pragmatic, as described below.

The issues in terms of biological mechanisms and computational benefits are elaborated below.

## Biological Mechanisms

Aside from the process of accurately capturing the nature of the biochemistry driving synaptic plasticity, there are some important questions that the biology may help us understand.

For example, a major unsolved problem is: how does the system know it is in the minus phase vs. the plus phase?

We start with the following assumptions based on existing work:

* The network dynamics have a natural theta-phase activity dynamic, via pulvinar predictive learning or hippocampal theta, which establishes a temporal ordering of minus-then-plus phases, at roughly the 200 msec theta cycle window.

* LTD / DAPK1 has a longer time constant of integration (medium scale in the XCAL framework), while LTP / CaMKII is faster (short scale in XCAL).  This aligns with the theta activation dynamics, such that DAPK1 naturally reflects more of the minus phase, while CaMKII reflects more of the plus phase.  In the current axon model, STau is 10 while MTau is 40, building on a common SSTau signal integrated with a 40msec time constant.

If you just learn based on the continuous kinase signals in a fully continuous manner, it probably would not result in reliable EDL signal (todo: test this!).  Thus, some additional "learn now" signal is likely required.  Here are some ideas.

* Learn at silence: due to sparse activity (10-20%), each synapse is relatively unlikely to experience continuous activity across both sender and receiver within a given theta cycle (e.g., $.1 * .1 = .01 = 1 / 100$).  Thus, there is likely to be a distinct window of synaptic activity, followed by silence.  This silence provides a reasonable timing cue for the end of the plus phase -- e.g., whatever happened roughly 200 msec prior to silence counts for learning.  The only problem would be if the plus phase goes completely silent, but this is likely to be relatively rare -- usually it is just a decrease in activity.  Also, there is a question about whether activity strictly on either side (pre only or post only) would cause a false alarm?  This relates to the Ca driver question below.

* The learn-at-silence idea also has the benefit of naturally blocking learning on repeated inputs -- only the last event in a sequence of repeated synaptic activations would lead to learning.

* The learn-at-silence idea is consistent with the Urakubo model, which targets the induction, not maintenance, functions of CaMKII: thus, there is some window after a bout of synaptic activity when these other actin-based functions of CaMKII kick in -- we need to account for that, and associated thresholds etc.   In particular, the system may be in a kind of labile kinase-battle state until this critical off signal occurs, at which point the relative balance between CaMKII and DAPK1 determines the direction of wt change. 

## Computational Issues: Average of Product vs. Product of Averages

One of the most important computational features of the Kinase learning rule is that it should be based on a local synaptic calcium signal that reflects the interaction of both pre and postsynaptic factors, integrated over time, in a way that could reflect detailed spike timing relationships between these two neurons.  By contrast, in Leabra and initially in Axon, the average activations were computed separately for sender and receiver, and then their product used for the CHL learning mechanism (i.e., Kinase should use an average of a product, vs. a product of averages in CHL).  In the rate-code based Leabra model, this did not matter too much because the neuron activations didn't change too much over the course of a trial, and there was no detailed spike timing information in the first place.

With discrete spiking, the *only* way to be sensitive to more detailed correlations in pre-post spiking activity is to compute the average of the products!  This entails significant computational cost in updating the average on a per-synapse, per-spike basis, but at least we should be able to do it per-spike and not every cycle, using the ISI values on Send and Recv neurons to optimize computation.

# KinaseAMax Abstract model

See https://github.com/emer/axon/tree/master/examples/kinaseq for a parallel simulation examining this abstract formulation and how to efficiently implement it.

The Leabra and initial Axon equations used cascading exponentially-updated running average variables on each neuron (recv and send) separately, with the final values multiplied for the CHL-like DWt function.  Here are the new names for these variables in the Kinase framework:

* `CaM` = first level integration: `CaM += (g * Spike - CaM) / MTau` -- spike drives up toward a max per-spike amount with gain factor (g = 8 std).  The g multiplier gives a faster effective rate constant for going up vs. decay back down -- this just sets the overall scale of the values.  This represents an abstract version of Ca, and CaM (calmodulin) which is activated by calcium, as driven by influx from NMDA and VGCC receptors.  This is `AvgSS` in Leabra.

* `CaP` = LTP, plus phase faster integration of CaM, reflecting CaMKII biologically.  Faster time constant makes this reflect the plus phase relative to the `CaD` LTD, minus-phase variable.  This is `AvgS` in Leabra.

* `CaD` = LTD, minus phase slower integration of CaP (cascaded one further step), reflecting DAPK1 biologically.  This is `AvgM` in Leabra.

The CHL product-based learning rule is a function of the difference between the plus-phase vs. minus-phase products of factors for S = Sender and R = Receiver:

* `DWt = SCaP * RCaP - SCaD * RCaD`

For the rate-code activations in Leabra, the product of these averages is likely to be similar to the average of the products at a synapse level, and computing neuron-level values is *much* faster computationally than integrating the products at the synapse level.  Indeed, experiments (a long time ago) showed no advantages to doing the synapse-level integration in Leabra.

# KinaseB*: Working up by simplifying the Urakubo model

To tackle the bottom-up approach toward deriving Kinase learning rules, we need to dramatically simplify the highly detailed Urakubo mechanisms into something that can be computed efficiently in large scale models.  This has turned out to be *much* simpler than originally estimated, with very simple exponential integration equations giving remarkably accurate approximations.  Furthermore, the critical synaptic Ca signal can be computed via an entirely factorized set of equations that depend separately on pre and post spiking, which is surprising.

## Allosteric NMDA and resulting Ca flux

![NMDAo in Urakubo vs. Kinase](results/fig_kinase_nmdao_thetaerr_nrep3_isi01.png?raw=true "Urakubo 08 NMDAo vs. simple Kinase NMDAo computation.")

First, the number of open NMDA channels is driven entirely by presynaptic spiking & Glu release, with an inhibitory factor, which can be very closely captured with simple exponential rate equations as shown in the above figure.  Given the significant complexity of the Urakubo allosteric NMDA model, it is remarkable how accurately a very simple model can capture it:

```Go
	ks.NMDAo -= ks.NMDAo / kp.DecayTau
	ks.NMDAi -= ks.NMDAi / kp.InhibTau
	if nex.PreSpike > 0 {
		ks.NMDAo += (1 - ks.NMDAi) * (kp.PreOpen - ks.NMDAo)
		ks.NMDAi += (1 - ks.NMDAi) * kp.PreInhib
	}
```

![Kinase Ca vs. Urakubo, VGCC0](results/fig_kinase_kinca_thetaerr_nrep3_isi01_vgcc0.png?raw=true "Kinase computed Ca vs. Urakubo 08, VGCC=0")

The presynaptically-driven NMDAo value can be combined with the standard voltage-gated dependency of the postsynaptic side (Mg blockage), to compute a net Ca value at each point in time, as shown in the above figure, which also fits the full Urakubo value (`PSD_Ca` = postsynaptic density Ca) very well.  The Vdrive factor is also important for this close fit -- this multiplies the Ca value by the same exponential function of voltage used in computing the VGCC conductance.

```Go
	ks.Ca = kp.CaGain * ks.NMDAo * nex.NMDAGmg
	if kp.CaVdrive {
		ks.Ca *= kp.Vdrive(chans.VToBio(vmd))
	}
	ks.Ca += kp.VGCCGain * nex.Gvgcc
```

![Kinase Ca vs. Urakubo, VGCC12](results/fig_kinase_kinca_thetaerr_nrep3_isi01_vgcc12.png?raw=true "Kinase computed Ca vs. Urakubo 08, VGCC=.12")

The previous figure shows the case where the VGCC currents are off, while the above one shows GBar for VGCC = .12, which is that replicates the Urakubo results.  A `VGCCCGain` of 20 provides a good fit to the overall resulting Ca spikes.

Thus, the bottom line is that a realistic synaptic Ca level can be computed *at any point* based on a presynaptically-driven `NMDAo` variable, and postsynaptically-driven `NMDAmg`, `Vdrive`, and `VGCC` factors.

## CaMKIIact from Ca

![Kinase CaP vs. Urakubo CaMKII](results/fig_kinase_cap_camkii_thetaerr_nrep1_isi08.png?raw=true "Kinase CaP vs. Urakubo 08")

The above figure shows that the response of CaMKIIact (activated CaMKII, in the PSD) can be accurately approximated via a simple cascaded exponential integration function of the form used in the existing Leabra and Axon models as described above.  Specifically, the CaP value shown in the plot is integrated with a `SpikeG = 90`, `MTau = 40`, and `PTau = 1200` (using `SynNMDACa` which just applies the cascaded integration dynamics to the raw `PSD_Ca` signal).

The `MTau = 40` value which effectively tracks the rise time of the curve is identical to that used in existing Axon models, except it is used for PTau, on top of an initial `MTau = 10` that roughly approximates the time dynamics of the Ca signal, based on raw input spiking.  Thus, overall, there is a good match for the rise time dynamics of this integration of CaMKII to the abstract model, which does effective error-driven learning at the 200 msec Theta timescale.

However, the decay time constant is dominated by the `PTau = 1200` and it is *much* longer than what is used in the model.  However, biologically, this likely reflects the fact that CaMKII persists in the N2B binding state to drive further steps of the plasticity process, and the most relevant factor here is the rise dynamics, and, really, the detailed relationship with DAPK1, as it is only the relative strength of these two that determines the direction of weight change.

In short, as with the earlier steps, it is remarkably straightforward to capture the dynamics of CaMKII integrating directly from the Ca influx signal, through the same cascaded exponential integration functions used in our effective computational models, even with the same relevant time dynamics.

All of this means that we can bridge directly to more abstract, computationally efficient formulations, which is important given the relative uncertainty about how DAPK1 actually works at the detailed biochemical level, and about how it interacts with CaMKII to drive DWt.  

Because we have concrete data from the Zito lab experiment confirming the theta window CHL-like learning dynamics, we know what the net result of the biochemistry is, and we have abstract models that capture it, so this greatly constraints from the "top-down" the exploration for how DAPK1 should behave.

## The VmDend wildcard

Despite all the positive developments above, there is a major unresolved issue in terms of understanding the actual behavior of the `Vm` membrane potential in the spine, which then significantly affects the postsynaptic NMDA and VGCC conductances.

TODO: explain and deal with this.  Issue is NMDA (and GABAB) for bistability vs. learning.

First step: just use Vm instead of VmDend (with a flag), and see how that does using full bio abstracted Ca signals per above.

## KinaseB1: SynNMDACa

The above approximations are all implemented efficiently through separate pre and postsynaptic variables, that are multiplied at the synaptic level to generate an accurate reflection of the synaptic Ca signal as would be given in the Urakubo model.  The `SynNMDACa` learning rule (aka KinaseB1) then simply runs the cascading integrators for CaP and CaD on top of this Ca signal, to get a learning rule.  Note that the first stage CaM signal is intended to capture Ca itself, so the time constant of integration for this level should be 1, when using this biologically-based Ca signal.

![SynNMDACa PSD_Ca thetaerr learning](results/fig_kinase_synnmdaca_thetaerr_nrep100_isi08_psdca.png?raw=true "SynNMDACa learning behavior, driven by PSD_Ca, in ThetaErrSweep @ 100 reps, SpikeG = 1, MTau = 1, PTau = 40, DTau = 40, DScale = 1.")

The above figure shows overall error-driven learning behavior of the SynNMDACa version, using the full Urakubo `PSD_Ca` calcium signal, on the ThetaErrSweep case with NRep = 100, ISISec = 0.8, VGCC Gbar = .12, SpikeG = 1, MTau = 1, PTau = 40, DTau = 40, DScale = 1, Lrate = 0.02.  It shows significant deviations from the error-driven ordering in the 50-25 and 25-50 cases, and the 100-50 is not very different from the no-error cases -- below we see that much of this is due to the VGCC dynamics.

![SynNMDACa KinCa thetaerr learning](results/fig_kinase_synnmdaca_thetaerr_nrep100_isi08_kinca.png?raw=true "SynNMDACa learning behavior, driven by Kinase Ca, in ThetaErrSweep @ 100 reps, SpikeG = 1, MTau = 1, PTau = 40, DTau = 40, DScale = 0.93.")

The above figure shows the comparable case using kinase-based Ca instead of full Urakubo, with VGCC = 0.12 -- note that it does differ significantly, in particular in the err LTD cases, suggesting that the Ca approximation, though good, needs a bit more tuning.

![SynNMDACa PSD_Ca VGCC = 0, thetaerr learning](results/fig_kinase_synnmdaca_thetaerr_nrep100_isi08_psdca_vgcc0.png?raw=true "SynNMDACa learning behavior, driven by PSD_Ca, with no VGCC, in ThetaErrSweep @ 100 reps, SpikeG = 1, MTau = 1, PTau = 40, DTau = 40, DScale = 1.")

Above figure shows same configuration but VGCC Gbar = 0, with DScale = 0.93 to balance weight changes.  It is a bit different, but not radically so.  Note that PTau = 10 (MTau = 1, DTau = 40) does not produce systematic error-driven learning pattern -- 10 for Ca and 40 for PTau seems better..

![SynNMDACa KinCa VGCC = 0 thetaerr learning](results/fig_kinase_synnmdaca_thetaerr_nrep100_isi08_kinca_vgcc0.png?raw=true "SynNMDACa learning behavior, driven by Kinase Ca, with no VGCC, in ThetaErrSweep @ 100 reps, SpikeG = 1, MTau = 1, PTau = 40, DTau = 40, DScale = 0.93.")

Finally, above figure shows KinaseCa version with VGCC=0 -- even worse.

TODO: Ca transients may be part of the issue?

TODO: CaM = 1, Vm instead of VmDend!  Definitely try VGCC!

TODO: redo thetaerr sweep results!

TODO: kinca VGCC rgclamp looks a lot better -- compare on raw Ca values.


# KinaseAMax behavior

## SynSpkCaOR

The following plots show the behavior of the most abstract `SynSpkCaOR` synaptic version of the Kinase algorithm, in capturing Ca values as represented by the CaM first-stage integration over discrete pre- post spikes.

This rule simply says that there is a "synaptic" spike impulse whenever *either* the pre or post: `SynSpk = SSpk || RSpk` -- either spike counts, but there is no specific interaction -- this is the least product-like.  This SynSpk value then drives the same cascade of time integrations, as usual, with the first CaM stage with a Tau of 10 serving to give a reasonable approximation of the biological Ca as computed above:

![NMDA / Ca vs. OR rule](results/fig_kinase_synspkca_or_thetaerr_nrep3_isi01.png?raw=true "Urakubo 08 PSD Ca vs. simple OR model of pre-post spike intergration in CaM signal.")

![NMDA / Ca vs. OR rule](results/fig_kinase_synspkca_or_thetaerr_nrep3_isi01_zoom.png?raw=true "Urakubo 08 PSD Ca vs. simple OR model of pre-post spike intergration in CaM signal.")

![NMDA / Ca vs. OR rule](results/fig_kinase_synspkca_or_thetaerr_nrep3_isi01_zoom2.png?raw=true "Urakubo 08 PSD Ca vs. simple OR model of pre-post spike intergration in CaM signal.")

![NMDA / Ca vs. OR rule](results/fig_kinase_synspkca_or_thetaerr_nrep3_isi01_zoom3.png?raw=true "Urakubo 08 PSD Ca vs. simple OR model of pre-post spike intergration in CaM signal.")

The above figures show that a very simple "OR" spike-driven Ca integration rule can *sort of* capture the complex allosteric dynamics from Urakubo.  The OR rule says that either a pre OR post spike drives an influx of Ca up to a max Ca driving level.  The CaM Tau = 10 here, with ThetaErr windows at 50 then 25 hz, repeated 3x with .1 sec ISI intervals.

TODO: point to thetaerr sweep results!

## RA25, Objrec, LVis sim results

This is the most abstract, pragmatic starting model, using same SS, S, M cascading running-average computation on top of synaptic Ca computed directly from NMDA channels and VmDend, etc (same as used for driving good activation dynamics).  Learning automatically happens at end of ThetaCycle.

Versions are noted by their Go emer/axon tag (bumping to 1.3.1 as start of the Kinase era).

## v1.3.1

* Urakubo-matching allosteric NMDA dynamics do NOT work well for soft bistability dynamic -- Inhib causes too much noise, like synaptic failure.  May work well at larger scale.  Sticking with `Tau = 100, MgC = 1, ITau = 1` (disabled).  Also `SnmdaDeplete` and `SeiDeplete` off -- also major source of noise, as discussed here: https://github.com/emer/axon/issues/28#issuecomment-1032297427

* First thing I tried worked really well: Basic default `40/10/40` taus for `CaM, CaP, CaD` cascade, `DScale=0.93`, `MinLrn=0.02`, `Lrate.Base=0.005`.

* `Jca` is pretty diffuse across all neurons -- not very discriminative from the recv side.  natural magnitudes around 5-10, which then translates to similar scale in Ca.

* `SnmdaO` is highly distinctive -- filter on this for winnowing synapses to update (very convenient given sender-based approach).  Need to init Ca's all to 0 in final WtFmDWt, exclude updating Ca for anything with straight below-thr on sender AvgS, AvgM.  Note that filtering based on ISI = -1 did not seem to work -- excludes first spike data.

* Overall, the pattern of weight change across recv neurons is very homogenous -- really need more of a postsynaptic selection factor.  VGCC?

* `lvis/sims/objrec` is showing significant hogging -- learns quickly then dies with hogging.  The lack of differentiation on recv side is clearly an issue.  

* Some issues: 

    + Not doing any actual time-integration of Ca -- all driven by instantaneous signals -- so VGCC is not really integrating that much over time -- only on CaM for 1 cycle.  In the biology, there is significant buffering / decay of Ca that would cause weaker signals to diminish significantly.  This is not captured.

    + `VmDend` is a major issue too -- stays relatively high, reducing differentiation, MgG is kind of high baseline level for everything.
    
    + About 10x slower overall performance-wise -- 
    
## v1.3.2

* Added `Dend.CaThr` on post-norm Ca (renamed Jca -> RCa) -- creates more differentiation on the recv side.  Reflects the effects of buffering.

* Optimization: 


# Urakubo

* Urakubo, H., Honda, M., Froemke, R. C., & Kuroda, S. (2008). Requirement of an allosteric kinetics of NMDA receptors for spike timing-dependent plasticity. The Journal of Neuroscience, 28(13), 3310–3323. http://www.ncbi.nlm.nih.gov/pubmed/18367598 | [Main PDF](https://github.com/emer/axon/blob/master/examples/urakubo/papers/UrakuboHondaFroemkeEtAl08.pdf) | [Supplemental PDF](https://github.com/emer/axon/blob/master/examples/urakubo/papers/UrakuboHondaFroemkeEtAl08_suppl.pdf) | [Model Details](https://github.com/emer/axon/blob/master/examples/urakubo/papers/UrakuboHondaFroemkeEtAl08_model_sup.pdf) (last one has all the key details but a few typos and the code is needed to get all the details right).

This model captures the complex postsynaptic chemical cascades that result in the changing levels of AMPA receptor numbers in the postsynaptic density (PSD) of a spine, as a function of the changing levels of calcium (Ca2+) entering the spine through NMDA and VGCC (voltage-gated calcium) channels.  The model mechanisms are strongly constrained bottom-up by the known kinetics of these chemical processes, and amazingly it captures many experimental results in the LTP / LTD experimental literature.  Thus, it provides an extremely valuable resource for understanding how plasticity works and how realistic but much simpler models can be constructed, which are actually usable in larger-scale simulations.

The complexity of this model of a single synaptic spine involves 116 different variables updated according to interacting differential equations requiring a time constant of 5e-5 (20 steps per the typical msec step used in axon).  In previous work with the original Genesis implementation of this model (see [CCN textbook](https://CompCogNeuro.org)), we were able to capture the net effect of this model on synaptic changes, at a rate-code level using random poisson spike trains, using the simple "checkmark" linearized version of the BCM learning rule (called "XCAL" for historical reasons), with an `r` value of around .9 (81% of the variance) -- see figure below.  However, with the actual spiking present in axon, there is an important remaining question as to whether a more realistic (yet still simplified) learning mechanism might be able to take advantage of some aspects of the detailed spiking patterns to do "better" somehow than the simplified rate code model.

The original Genesis model is a good object lesson in how a generic modeling framework can obscure what is actually being computed.  The model shows up as a complex soup of crossing lines in the Genesis GUI, and the code is a massive dump of inscrutable numbers, with key bits of actual computation carefully hidden across many different bits of code, with seemingly endless flags and conditions making it very difficult to determine what code path is actually in effect.  In addition, Genesis is ancient software that requires an old gcc compiler to build it, and it only runs under X windows.

By contrast, the present model is purpose-built in compact, well-documented highly readable Go code, using a few very simple types in the [emergent chem](https://github.com/emer/emergent/tree/master/chem) package, providing the basic `React` and `Enz` mechanisms used for simulating the chemical reactions.

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


# DAPK1 and Competitive NMDA GluN2B binding of CaMKII for LTD

There is recent evidence for a competitive interaction between DAPK1 and CaMKII binding at the GluN2B subunit of the NMDA receptor (Goodell et al., 2017), which is important for the LFS LTD, that the model does not currently capture.  Thus, it is important to add this dynamic, tuning it to capture the LFS LTD effects, and then testing whether this is sufficient to capture the CHL error-driven learning dynamic of primary interest.

## CaMKII binding to GluN2B

![CaMKII GluNR2B Binding](figs/fig_urakubo_camkii_nums.png?raw=true "CaMKII Ca/CaM binding and autophosphorylation, showing the NR2B binding modulation at reactions 6 and 9, in terms of reduced backward rate constant 0.02 instead of 400, 0.001 and 1, along with 8 going forward at 9.  Adapted from SI3 in Urakubo modeling suppl.")

The Urakubo model does not directly include the role of CaMKII binding to the GluN2B subunit of the NMDAR (Bayer et al, 2001; Coultrap & Bayer, 2012; Barcomb et al, 2016), but does have indirect accounting for it.  There are two primary effects:

* Concentrating CaMKII at the PSD on the NMDA receptors depends on GluN2B binding.  Only the PSD-localized CaMKII can activate and stabilize AMPAR in the PSD, so this is the primary determinant of LTP effects of CaMKII.  One site of GluN2B binding is dependent on the CaMKII T286 autophosphorylation (CaMKIIP in the model), which is in turn a complex function of Ca-CaM binding, driven by Ca influx, and another is dependent on Ca-CaM binding (Bayer et al, 2001).  This is reflected in the model by the asymmetric diffusion constants for CaMKIIP, which is 10 times larger for Cyt -> PSD than the PSD -> Cyt direction (6 vs. 0.6).

* Altering the CaM binding constants: "modified the binding constants between bare CaM and the inactive subunit in consideration of experimentally observed basal activity (Fukunaga et al., 1993; Kawaguchi and Hirano, 2002). This modification also implicitly represents CaMKII by binding to NMDARs (Bayer et al., 2001; Leonard et al., 2002)".  In the adapted Figure SI3 shown above, it is reactions 4 and 9: 4 is CaM + plain CaMKII (0.0004 Kf, 1 Kb), 9 is CaM + CaMKIIP (8 Kf, 0.001 Kb).

Thus, to incorporate DAPK1, we first need to simulate the GluN2B binding dynamics explicitly.

Shen & Meyer (1999), Bayer et al (2001) and Leonard et al (2002) have the relevant details on GluN2B binding.  The key points are:

* There are two N2B binding sites: a -C and -P, and the -C is determined by Ca/CaM activation of CaMKII (only), while -P depends on T286 autophosphorylation.  In general, the dynamics are that Ca/CaM activation causes the initial binding, and T286 keeps it around longer after Ca/CaM is gone. 

* >In fact, this NR2B domain contains two sites with different modes of regulated CaMKII binding (Fig. 1), a Ca2+/CaM-regulated site within residues 1,120 -- 1,482 of NR2B (NR2B-C, later determined to be 1303) and a phosphorylation-regulated site within residues 839 -- 1,120 (NR2B-P). 

* >During washes with Ca2+ but without CaM, however, CaM remained bound to the complex. Thus, interaction with NR2B-C seems to increase the CaM affinity of CaMKII, similar to CaM trapping by the autophosphorylated kinase, as CaM dissociates from unphosphorylated kinase that is not complexed with NR2B-C with an off rate of approximately 2 s-1 (ref. 11)

* >Immunoblotting with an antibody sensitive to phosphorylated T286 showed that the autonomously active NR2B-C-bound enzyme is not phosphorylated at T286.

Also, the Urakubo model centrally features Ca-CaM binding to NMDA receptors to mediate the allosteric effects -- this is not on GluN2B, but rather on the NR1 C0 region -- which is on the 2A subunit, not 2B -- so an entirely different mechanism (e.g., Akyol et al, 04).  Interestingly, the NR2A subtype is only expressed later in development, so the allosteric effect may have a developmental timecourse (Umemiya et al, 2001).

In preparation for understanding DAPK1, here's a summary of the key AutoP dynamic on CaMKII, as summarized in above figure.  First, T286 P and Ca/CaM binding are *separate processes* that nevertheless interact strongly.  There are separate vars for for P, noP, and various levels of Ca/CaM binding.  Further there is a "total activity" factor that integrates across P and Ca/CaM states, and yet another GluN2B binding function, so there are 4 primary dynamics:

* **CaM** -- CaM is attracted to T286P and binds strongly to bare CaMKIIP (eq 9, Kf=8, Kb=0.001), but otherwise it takes Ca+ activation to bind to CaMKII (noP) (eq 5, Kf=8 -- only for N2B bound, otherwise inactive binding rates).

only binds to noP CaMKII, *except* for the very first binding to CaMKIIP (eq 9), which has a high Kf (8) and slow Kb (0.001).  P *or* N2B is responsible for this high CaM binding, so eq 9 should apply to either case.  Any N2B bound CaMKII (noP) should have same as 9 affinity.  Also, T286P prevents further CaM binding (presumably the CaM is actually bound to T286?), and P also reduces magnitude of Ca binding to CaM (eq 8) but balances the Kf, Kb (both 1).  Question: should N2B also induce same 8 dynamics in CaMKII plain case?

* **Auto.K** -- determines rate of further AutoP, and it is a complicated function of current P and Ca/CaM binding.  Critically, noP but Ca/CaM binding contributes positively to rate of AutoP -- not only P itself.  And there is some "mass action" effect that is a function of any Ca/CaM or P state.  More P does accelerate this function, but not massively so: if completely left alone and driven by constant Ca input, it accelerates, but PP1 quickly turns off.  Also Ca/CaM unbinds with decent speed.

* **Auto.Act** -- total activation: to first order, any CaMKII with *either* P or Ca/CaM binding contributes, with weak diffs (e.g., 3CaCaMP is the most active, contributes with a factor of 1, but other ones contribute with factors between .75 - .8, so not that big of a diff).  This total act is what drives the further kinase activity for AMPAR trafficking.

* **GluN2B** binding -- There is a *distinct* function for GluN2B binding.  Per above, `GluN2BCaCaM` (10 Kf, 0.001 Kb) drives N2B binding as function of any CaM (does not accelerate for simplicity).  `GluN2BP` for CaMKIIP, no CaM (0.01 Kf, 0.0001 Kb) keeps P bound (low Kb), but does not drive much forward binding.  `GluN2NoP` for CaMKII (noP) (0.00001 Kf, 1000 Kb) unbinds noP (high Kb) quickly, and does not do much forward binding (todo just set to 0?)

## DAPK1 dynamics

![DAPK1 GluNR2B Binding](figs/fig_urakubo_dapk1.png?raw=true "DAPK1 Ca/CaM binding and autophosphorylation, showing the NR2B binding modulation at reactions 6 and 9, in terms of reduced backward rate constant 0.02 instead of 400, 0.001 and 1, along with 8 going forward at 9.  Adapted from SI3 in Urakubo modeling suppl.")

Goodell et al. (2017) provide a detailed investigation of DAPK1 dynamics.

* DAPK1/GluN2B binding is homologous to the well-studied CaMKII/GluN2B binding to the same region.

* However, it has the opposite dependence on Ca-CaM binding: inversely proportional -- key question is exactly how.

* DAPK1 autophosphorylates at the Ser308 site (instead of Thr286, which is CaMKII spot).  de Diego et al., (2010) suggest that the S308 site is similar to T286 in having multiple additional binding steps, and remains phosphorylated after CaM dissociation (autophosphorylation).  Thus, it is reasonable to apply the same complicated T286 logic here it seems, but there are some twists.

* Per Shani et al. (2001), "It was found that in the absence of Ca2+/CaM, DRP-1 displayed unexpectedly high basal levels of autophosphorylation. The latter was strongly inhibited by Ca2+/CaM concentrations that fully activate the enzyme, suggesting inverse relationships between autophosphorylation and substrate phosphorylation."

    DAPK1P (S380 AutoP) = inactive form for binding to GluN2B (akin to its substrate P. effects) -- this is default state.  Data: mutant S308A cannot autoP, and shows increased substrate binding / apoptosis.  mutant S308D is always autoP (aspartic acid), and prevented apoptosis.

* Per Shohat et al (2001), "The autophosphorylation of Ser308 was found to be Ca2/CaM-independent; moreover, it was strongly inhibited by the addition of Ca2/CaM."  Thus, Ca/CaM inhibits P, and P inhibits Ca/CaM reciprocally.

* In summary, it starts out fully P, then with some initial deP from CaN, further Ca/CaM binding further inhibits additional autoP pressure -- so it integrates Ca/CaM into levels of deP and its overall kinetic binding activity is a function of Ca/CaM + deP (this can be a mirror of CaMKII where it is Ca/CaM + P).  Assuming strong binding is due to active Ca3CaM form of CaM -- no indication that raw CaM is bound?

* There does not appear to be a strong interaction between CaM binding and GluN2B binding?  Most of the interaction is with AutoP causing P, and CaM inhibiting that.  But GluN2B binding is just more-or-less an outcome of these CaM, deP vs. AutoP dynamics, without changing them significantly?  TODO: verify this!

* Goodell et al (2017) show that deP is due to CaN, although other lit shows PP2A deP activity, but that is likely in non-neuronal context (Bialik & Kimchi, 2014), and PP2A is not in PSD.

* To implement, we need Act to have reversed roles of P and noP, but AutoK is still driven by P but inhibited by Ca/CaM levels (invert that for example).

* Thus, overall, DAPK1 is default baseline autophosphorylated.  And it is relatively insensitive to Ca/CaM in this state, because of the inhibited binding constants.  However, when you get some initial CaN dephosphorylation activity, that then allows more Ca-CaM binding, which further inhibits autoP, and you move further away from full autoP.  It is this noP DAPK1 that is actually "active" catalytically: this is when it binds to GluN2B most strongly.

* So, key point: CaN is "pulling the trigger" for a learning event, after which time DAPK1 effectively measures / Camintegrates the Ca/CaM through a kind of inverse-autophosphorylation process (the opposite form of positive feedback loop present in CaMKII).  This is then the "opponent process" to the CaMKII, both of which are competing for GluN2B.  If DAPK1 wins out, then LTD happens.  If CaMKII wins, then LTP happens.  But critically, both are directly integrating Ca signals through a similar integrative, rate-sensitive dynamic!  However, *if* DAPK1 preferentially reflects the earlier minus phase time frame of activity (i.e., it decays more quickly), then it directly represents this minus phase co-product, while CaMKII gets to be the plus phase.  It could really be that simple and direct!

* CaN is the primary source of dephosphorylation at S308, such that CaN is primary initiator of DAPK1 binding to GluN2B.  CaN is directly activated by Ca/CaM, so it is in a good position to detect the onset of a learning window, opening up the DAPK1 to be sensitive to further Ca levels.  It is unclear (yet) if PP1 and PP2A, which dephosphorylate the T286 site on CaMKII, might also play a similar role on DAPK1, and if not, why not?  PP1 does de-activate CaN however, so indirectly it is having that effect, and it may thus set the overall time window for a given plasticity event.

* In summary, at the most abstract, computational level, CaMKII represents the plus phase $<xy>_s$ signal, and DAPK1 represents the minus phase $<xy>_m$ signal, in terms of time-integrated Ca++ signals across these phases. These two kinases share many properties, but differ critically in their signs (at multiple levels), and this makes them work just like a classic opponent-process: Go - NoGo, excitation vs. inhibition, etc. They both compete at the NMDA GluN2B site, so that is the final locus of the opponent competitive dynamic. If DAPK1 wins, LTD happens. If CaMKII wins, LTP happens. If there is a tie, nothing happens! Perfect error-driven learning! Furthermore, both rely on a similar auto-phosphorylation (autoP) dynamic to integrate over recent synaptic activity.

# TODO

* restrict GluN2BN?

* more explicit division relationship between CaMKII and DAPK1?

* add an abstract model with just integration of Ca at 2 diff rate constants, as in current XCAL, and compare that -- time to try a more top-down approach again to see what that looks like. To first approximation, DAPK1 and CaMKII are just integrating Ca with different dynamics..

# References

* Akyol, Z., Bartos, J. A., Merrill, M. A., Faga, L. A., Jaren, O. R., Shea, M. A., & Hell, J. W. (2004). Apo-Calmodulin Binds with its C-terminal Domain to the N-Methyl-d-aspartate Receptor NR1 C0 Region *. Journal of Biological Chemistry, 279(3), 2166–2175. https://doi.org/10.1074/jbc.M302542200

* Barcomb, K., Hell, J. W., Benke, T. A., & Bayer, K. U. (2016). The CaMKII/GluN2B Protein Interaction Maintains Synaptic Strength. *Journal of Biological Chemistry, 291(31),* 16082–16089. https://doi.org/10.1074/jbc.M116.734822

* Bayer, K.-U., De Koninck, P., Leonard, A. S., Hell, J. W., & Schulman, H. (2001). Interaction with the NMDA receptor locks CaMKII in an active conformation. *Nature, 411(6839),* 801–805. https://doi.org/10.1038/35081080

* Bhalla, U. S., & Iyengar, R. (1999). Emergent properties of networks of biological signaling pathways. *Science.* https://doi.org/10.1126/science.283.5400.381

* Coultrap, S. J., & Bayer, K. U. (2012). CaMKII regulation in information processing and storage. *Trends in Neurosciences, 35(10),* 607–618. https://doi.org/10.1016/j.tins.2012.05.003

* de Diego, I., Kuper, J., Bakalova, N., Kursula, P., & Wilmanns, M. (2010). Molecular basis of the death-associated protein kinase-calcium/calmodulin regulator complex. Science Signaling, 3(106), ra6. https://doi.org/10.1126/scisignal.2000552

* Dudek, S. M., & Bear, M. F. (1992). Homosynaptic long-term depression in area CA1 of hippocampus and effects of N-methyl-D-aspartate receptor blockade. *Proceedings of the National Academy of Sciences U. S. A., 89(10),* 4363–4367. http://www.ncbi.nlm.nih.gov/pubmed/1350090

* Doi, T., Kuroda, S., Michikawa, T., & Kawato, M. (2005). Inositol 1,4,5-trisphosphate-dependent Ca2+ threshold dynamics detect spike timing in cerebellar purkinje cells. *Journal of Neuroscience, 25(4),* 950–961. https://doi.org/10.1523/JNEUROSCI.2727-04.2005

* Dupont, G., Houart, G., & De Koninck, P. (2003). Sensitivity of CaM kinase II to the frequency of Ca2+ oscillations: A simple model. *Cell Calcium, 34(6),* 485–497. https://doi.org/10.1016/S0143-4160(03)00152-0

* Goodell, D. J., Zaegel, V., Coultrap, S. J., Hell, J. W., & Bayer, K. U. (2017). DAPK1 mediates LTD by making CaMKII/GluN2B binding LTP specific. *Cell Reports, 19(11),* 2231–2243. https://doi.org/10.1016/j.celrep.2017.05.068

* Kuroda, S., Schweighofer, N., & Kawato, M. (2001). Exploration of signal transduction pathways in cerebellar long-term depression by kinetic simulation. *Journal of Neuroscience, 21(15),* 5693–5702. https://doi.org/10.1523/JNEUROSCI.21-15-05693.2001

* Migliore, M., Hoffman, D. A., Magee, J. C., & Johnston, D. (1999). Role of an A-Type K+ Conductance in the Back-Propagation of Action Potentials in the Dendrites of Hippocampal Pyramidal Neurons. *Journal of Computational Neuroscience, 7(1),* 5–15. https://doi.org/10.1023/A:1008906225285

* Poirazi, P., Brannon, T., & Mel, B. W. (2003). Arithmetic of Subthreshold Synaptic Summation in a Model CA1 Pyramidal Cell. *Neuron, 37(6),* 977–987. https://doi.org/10.1016/S0896-6273(03)00148-X

* Shani, G., Henis-Korenblit, S., Jona, G., Gileadi, O., Eisenstein, M., Ziv, T., Admon, A., & Kimchi, A. (2001). Autophosphorylation restrains the apoptotic activity of DRP-1 kinase by controlling dimerization and calmodulin binding. The EMBO Journal, 20(5), 1099–1113. https://doi.org/10.1093/emboj/20.5.1099

* Urakubo, H., Honda, M., Froemke, R. C., & Kuroda, S. (2008). Requirement of an allosteric kinetics of NMDA receptors for spike timing-dependent plasticity. *The Journal of Neuroscience, 28(13),* 3310–3323. http://www.ncbi.nlm.nih.gov/pubmed/18367598

* Umemiya, M., Chen, N., Raymond, L. A., & Murphy, T. H. (2001). A Calcium-Dependent Feedback Mechanism Participates in Shaping Single NMDA Miniature EPSCs. Journal of Neuroscience, 21(1), 1–9. https://doi.org/10.1523/JNEUROSCI.21-01-00001.2001

