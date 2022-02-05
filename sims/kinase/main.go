// Copyright (c) 2021, The Emergent Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*
kinase: This simulation replicates the Kinase et al, 2008 detailed model of spike-driven
learning, including intracellular Ca-driven signaling, involving CaMKII, CaN, PKA, PP1.
*/
package main

import (
	"fmt"
	"log"
	"strconv"
	"strings"

	"github.com/emer/axon/axon"
	"github.com/emer/axon/chans"
	"github.com/emer/emergent/chem"
	"github.com/emer/emergent/emer"
	"github.com/emer/emergent/netview"
	"github.com/emer/emergent/params"
	"github.com/emer/etable/eplot"
	"github.com/emer/etable/etable"
	"github.com/emer/etable/etensor"
	_ "github.com/emer/etable/etview" // include to get gui views
	"github.com/goki/gi/gi"
	"github.com/goki/gi/gimain"
	"github.com/goki/gi/giv"
	"github.com/goki/ki/ki"
	"github.com/goki/ki/kit"
	"github.com/goki/mat32"
)

// this is the stub main for gogi that calls our actual mainrun function, at end of file
func main() {
	gimain.Main(func() {
		mainrun()
	})
}

// LogPrec is precision for saving float values in logs -- requires 6 not 4!
const LogPrec = 6

// SimOpts has high-level simulation options that are accessed in the code
type SimOpts struct {
	InitBaseline bool `def:"true" desc:"use 500 sec pre-compiled baseline for initialization"`
	UseN2B       bool `def:"true" desc:"use the GluN2B binding for CaMKII dynamics -- explicitly breaks out this binding and its consequences for localizing CaMKII in the PSD, but without UseDAPK1, it should replicate original Kinase dynamics, as it does not include any competition there."`
	UseDAPK1     bool `desc:"use the DAPK1 competitive GluN2B binding (departs from standard Kinase -- otherwise the same."`
}

// TheOpts are the global sim options
var TheOpts SimOpts

func (so *SimOpts) Defaults() {
	so.InitBaseline = true
	so.UseN2B = true
	so.UseDAPK1 = true
}

// ParamSets for basic parameters
// Base is always applied, and others can be optionally selected to apply on top of that
var ParamSets = params.Sets{
	{Name: "Base", Desc: "these are the best params", Sheets: params.Sheets{
		"Network": &params.Sheet{
			{Sel: "Layer", Desc: "all defaults",
				Params: params.Params{
					"Layer.Act.Spike.Tr":     "7",
					"Layer.Act.Spike.RTau":   "3", // maybe could go a bit wider even
					"Layer.Act.Dt.VmTau":     "1",
					"Layer.Act.Dt.VmDendTau": "1",
					"Layer.Act.Dt.VmSteps":   "2",
					"Layer.Act.Dt.GeTau":     "1",    // not natural but fits spike current injection
					"Layer.Act.VmRange.Max":  "0.97", // max for dendrite
					"Layer.Act.Spike.ExpThr": "0.9",  // note: critical to keep < Max!
					// Erev = .35 = -65 instead of -70
					"Layer.Act.Spike.Thr": ".55", // also bump up
					"Layer.Act.Spike.VmR": ".45",
					"Layer.Act.Init.Vm":   ".35",
					"Layer.Act.Erev.L":    ".35",
				}},
		},
	}},
}

// Extra state for neuron -- VGCC and AK
type NeuronEx struct {
	Gvgcc      float32 `desc:"VGCC total conductance"`
	VGCCm      float32 `desc:"VGCC M gate -- activates with increasing Vm"`
	VGCCh      float32 `desc:"VGCC H gate -- deactivates with increasing Vm"`
	VGCCJcaPSD float32 `desc:"VGCC Ca calcium contribution to PSD"`
	VGCCJcaCyt float32 `desc:"VGCC Ca calcium contribution to Cyt"`
	Gak        float32 `desc:"AK total conductance"`
	AKm        float32 `desc:"AK M gate -- activates with increasing Vm"`
	AKh        float32 `desc:"AK H gate -- deactivates with increasing Vm"`
	PreSpike   float32 `desc:"1 = the presynaptic neuron spiked"`
	PreSpikeT  float32 `desc:"time when pre last spiked, in sec (from spine.Time)"`
	PreISI     float32 `desc:"ISI between last spike and prior one"`
}

func (nex *NeuronEx) Init() {
	nex.Gvgcc = 0
	nex.VGCCm = 0
	nex.VGCCh = 1
	nex.VGCCJcaPSD = 0
	nex.VGCCJcaCyt = 0
	nex.Gak = 0
	nex.AKm = 0
	nex.AKh = 1
	nex.PreSpike = 0
	nex.PreSpikeT = -1
	nex.PreISI = -1
}

// Sim encapsulates the entire simulation model, and we define all the
// functionality as methods on this struct.  This structure keeps all relevant
// state information organized and available without having to pass everything around
// as arguments to methods, and provides the core GUI interface (note the view tags
// for the fields which provide hints to how things should be displayed).
type Sim struct {
	Net          *axon.Network            `view:"no-inline" desc:"the network -- click to view / edit parameters for layers, prjns, etc"`
	Spine        Spine                    `desc:"the spine state with Kinase intracellular model"`
	Neuron       *axon.Neuron             `view:"no-inline" desc:"the neuron"`
	NeuronEx     NeuronEx                 `view:"no-inline" desc:"extra neuron state for additional channels: VGCC, AK"`
	KinaseSyn    KinaseState              `view:"no-inline" desc:"kinase synapse state values"`
	KinaseParams KinaseParams             `view:"no-inline" desc:"kinase parameters"`
	Params       params.Sets              `view:"no-inline" desc:"full collection of param sets"`
	Stim         Stims                    `desc:"what stimulation to drive with"`
	KinaseOnly   bool                     `desc:"only run the kinase algorithm, not the detailed biophysical Urakobo-based model"`
	ISISec       float64                  `desc:"inter-stimulus-interval in seconds -- between reps"`
	NReps        int                      `desc:"number of repetitions -- takes 100 to produce classic STDP"`
	FinalSecs    float64                  `def:"20,50,100" desc:"number of seconds to run after the manipulation -- results are strongest after 100, decaying somewhat after that point -- 20 shows similar qualitative results but weaker, 50 is pretty close to 100 -- less than 20 not recommended."`
	DurMsec      int                      `desc:"duration for activity window"`
	SendHz       float32                  `desc:"sending firing frequency (used as minus phase for ThetaErr)"`
	RecvHz       float32                  `desc:"receiving firing frequency (used as plus phase for ThetaErr)"`
	GeStim       float32                  `desc:"stimulating current injection"`
	DeltaT       int                      `desc:"in msec, difference of Tpost - Tpre == pos = LTP, neg = LTD STDP"`
	DeltaTRange  int                      `desc:"range for sweep of DeltaT -- actual range is - to +"`
	DeltaTInc    int                      `desc:"increment for sweep of DeltaT"`
	RGClamp      bool                     `desc:"use Ge current clamping instead of distrete pulsing for firing rate-based manips, e.g., ThetaErr"`
	Opts         SimOpts                  `view:"inline" desc:"global simulation options controlling major differences in behavior"`
	GluN2BN      float64                  `desc:"total initial amount of GluN2B"`
	DAPK1AutoK   float64                  `desc:"strength of AutoK autophosphorylation of DAPK1 -- must be strong enough to balance CaM drive"`
	DAPK1_AMPAR  float64                  `desc:"strength of AMPAR inhibitory effect from DAPK1"`
	DAPK1lrate   float64                  `desc:"multiplier for diff between DAPK1 and CaMKII"`
	CaNDAPK1     float64                  `desc:"Km for the CaM dephosphorylation of DAPK1"`
	VmDend       bool                     `desc:"use dendritic Vm signal for driving spine channels"`
	NMDAAxon     bool                     `desc:"use the Axon NMDA channel instead of the allosteric Kinase one"`
	NMDAGbar     float32                  `def:"0,0.15" desc:"strength of NMDA current -- 0.15 default for posterior cortex"`
	GABABGbar    float32                  `def:"0,0.2" desc:"strength of GABAB current -- 0.2 default for posterior cortex"`
	VGCC         chans.VGCCParams         `desc:"VGCC parameters: set Gbar > 0 to include"`
	AK           chans.AKParams           `desc:"A-type potassium channel parameters: set Gbar > 0 to include"`
	CaTarg       CaState                  `desc:"target calcium level for CaTarg stim"`
	InitWt       float64                  `inactive:"+" desc:"initial weight value: Trp_AMPA value at baseline"`
	Logs         map[string]*etable.Table `view:"no-inline" desc:"all logs"`
	Plots        map[string]*eplot.Plot2D `view:"-" desc:"all plots"`

	// internal state - view:"-"
	Msec      int              `inactive:"+" desc:"current cycle of updating"`
	Win       *gi.Window       `view:"-" desc:"main GUI window"`
	NetView   *netview.NetView `view:"-" desc:"the network viewer"`
	ToolBar   *gi.ToolBar      `view:"-" desc:"the master toolbar"`
	IsRunning bool             `view:"-" desc:"true if sim is running"`
	StopNow   bool             `view:"-" desc:"flag to stop running"`
}

// this registers this Sim Type and gives it properties that e.g.,
// prompt for filename for save methods.
var KiT_Sim = kit.Types.AddType(&Sim{}, SimProps)

// TheSim is the overall state for this simulation
var TheSim Sim

// New creates new blank elements and initializes defaults
func (ss *Sim) New() {
	ss.Opts.Defaults()
	ss.Spine.Defaults()
	ss.Spine.Init()
	ss.KinaseSyn.Init()
	ss.InitWt = ss.Spine.States.AMPAR.Trp.Tot
	ss.Net = &axon.Network{}
	ss.Params = ParamSets
	ss.Stim = STDP // ThetaErrComp
	ss.ISISec = 0.8
	ss.NReps = 1
	ss.FinalSecs = 0
	ss.DurMsec = 200
	ss.SendHz = 50
	ss.RecvHz = 25
	ss.DeltaT = 16
	ss.DeltaTRange = 50
	ss.DeltaTInc = 5
	ss.RGClamp = true
	ss.Defaults()
}

// Defaults sets default params
func (ss *Sim) Defaults() {
	ss.Opts.Defaults()
	ss.Spine.Defaults()
	ss.KinaseParams.Defaults()
	ss.GluN2BN = 2
	ss.DAPK1AutoK = DAPK1AutoK
	ss.DAPK1_AMPAR = ss.Spine.AMPAR.Phos.DAPK1_AMPAR
	ss.DAPK1lrate = ss.Spine.AMPAR.Phos.DAPK1lrate
	ss.CaNDAPK1 = 11
	ss.GeStim = 2
	ss.NMDAGbar = 0.15 // 0.1 to 0.15 matches pre-spike increase in vm
	ss.GABABGbar = 0.0 // 0.2
	ss.VGCC.Defaults()
	ss.VGCC.Gbar = 0.12 // 0.12 matches vgcc jca
	ss.AK.Defaults()
	ss.AK.Gbar = 0 // todo: figure this out!
	ss.CaTarg.Cyt = 10
	ss.CaTarg.PSD = 10
}

func (ss *Sim) Log(name string) *etable.Table {
	if ss.Logs == nil {
		ss.Logs = make(map[string]*etable.Table)
	}
	dt, ok := ss.Logs[name]
	if ok {
		return dt
	}
	dt = &etable.Table{}
	ss.Logs[name] = dt
	return dt
}

func (ss *Sim) Plot(name string) *eplot.Plot2D {
	return ss.Plots[name]
}

func (ss *Sim) AddPlot(name string, plt *eplot.Plot2D) {
	if ss.Plots == nil {
		ss.Plots = make(map[string]*eplot.Plot2D)
	}
	ss.Plots[name] = plt
}

////////////////////////////////////////////////////////////////////////////////////////////
// 		Configs

// Config configures all the elements using the standard functions
func (ss *Sim) Config() {
	ss.ConfigNet(ss.Net)
	ss.ConfigDWtLog(ss.Log("DWtLog"))
	ss.ConfigPhaseDWtLog(ss.Log("PhaseDWtLog"))
	ss.ConfigTimeLog(ss.Log("MsecLog"))
	ss.ConfigTimeLog(ss.Log("Msec10Log"))
	ss.ConfigTimeLog(ss.Log("Msec100Log"))
	ss.ConfigTimeLog(ss.Log("MsecLog2"))
	ss.ConfigTimeLog(ss.Log("Msec10Log2"))
	ss.ConfigTimeLog(ss.Log("Msec100Log2"))
	ss.ConfigAutoKLog(ss.Log("AutoKLog"))
}

func (ss *Sim) ConfigNet(net *axon.Network) {
	net.InitName(net, "Neuron")
	ly := net.AddLayer2D("Neuron", 1, 1, emer.Hidden).(*axon.Layer)
	net.Defaults()
	err := net.Build()
	if err != nil {
		log.Println(err)
		return
	}
	ss.InitWts(net)
	ss.Neuron = &ly.Neurons[0]
}

// InitWts loads the saved weights
func (ss *Sim) InitWts(net *axon.Network) {
	net.InitWts()
}

// SetParams sets the params for "Base" and then current ParamSet.
// If sheet is empty, then it applies all avail sheets (e.g., Network, Sim)
// otherwise just the named sheet
// if setMsg = true then we output a message for each param that was set.
func (ss *Sim) SetParams(sheet string, setMsg bool) error {
	if sheet == "" {
		// this is important for catching typos and ensuring that all sheets can be used
		ss.Params.ValidateSheets([]string{"Network", "Sim"})
	}
	err := ss.SetParamsSet("Base", sheet, setMsg)
	// if ss.ParamSet != "" && ss.ParamSet != "Base" {
	// 	sps := strings.Fields(ss.ParamSet)
	// 	for _, ps := range sps {
	// 		err = ss.SetParamsSet(ps, sheet, setMsg)
	// 	}
	// }
	return err
}

// SetParamsSet sets the params for given params.Set name.
// If sheet is empty, then it applies all avail sheets (e.g., Network, Sim)
// otherwise just the named sheet
// if setMsg = true then we output a message for each param that was set.
func (ss *Sim) SetParamsSet(setNm string, sheet string, setMsg bool) error {
	pset, err := ss.Params.SetByNameTry(setNm)
	if err != nil {
		return err
	}
	if sheet == "" || sheet == "Network" {
		netp, ok := pset.Sheets["Network"]
		if ok {
			ss.Net.ApplyParams(netp, setMsg)
		}
	}

	if sheet == "" || sheet == "Sim" {
		simp, ok := pset.Sheets["Sim"]
		if ok {
			simp.Apply(ss, setMsg)
		}
	}
	// note: if you have more complex environments with parameters, definitely add
	// sheets for them, e.g., "TrainEnv", "TestEnv" etc
	return err
}

////////////////////////////////////////////////////////////////////////////////
// 	    Init, utils

// Init restarts the run, and initializes everything, including network weights
// and resets the epoch log table
func (ss *Sim) Init() {
	TheOpts = ss.Opts
	ss.Spine.AMPAR.Phos.DAPK1_AMPAR = ss.DAPK1_AMPAR
	ss.Spine.AMPAR.Phos.DAPK1lrate = ss.DAPK1lrate
	DAPK1AutoK = ss.DAPK1AutoK
	ss.Spine.NMDAR.GluN2BN = chem.CoToN(ss.GluN2BN, PSDVol)             // 120 works for N2B only (no DAPK1)
	ss.Spine.DAPK1.CaNSer308.SetKmVol(ss.CaNDAPK1, CytVol, 1.34, 0.335) // 10: 11 μM Km = 0.0031724
	ss.Spine.Init()
	ss.KinaseSyn.Init()
	ss.InitWt = ss.Spine.States.AMPAR.Trp.Tot
	ss.NeuronEx.Init()
	ss.Msec = 0
	ss.SetParams("", false) // all sheets
	ly := ss.Net.LayerByName("Neuron").(axon.AxonLayer).AsAxon()
	if ss.NMDAAxon {
		ly.Act.NMDA.Gbar = ss.NMDAGbar
	} else {
		ly.Act.NMDA.Gbar = 0
	}
	ly.Act.GABAB.Gbar = ss.GABABGbar
	ss.InitWts(ss.Net)
	ss.StopNow = false
	ss.UpdateView()
}

// Counters returns a string of the current counter state
// use tabs to achieve a reasonable formatting overall
// and add a few tabs at the end to allow for expansion..
func (ss *Sim) Counters() string {
	return fmt.Sprintf("Msec:\t%d\t\t\t", ss.Msec)
}

func (ss *Sim) UpdateView() {
	if ss.NetView != nil && ss.NetView.IsVisible() {
		ss.NetView.Record(ss.Counters())
		// note: essential to use Go version of update when called from another goroutine
		ss.NetView.GoUpdate() // note: using counters is significantly slower..
	}
}

////////////////////////////////////////////////////////////////////////////////
// 	    Running the Network, starting bottom-up..

// RunStim runs current Stim selection
func (ss *Sim) RunStim() {
	fn, has := StimFuncs[ss.Stim]
	if !has {
		fmt.Printf("Stim function: %s not found!\n", ss.Stim)
		return
	}
	ss.StopNow = false
	go fn()
}

// NeuronUpdt updates the neuron and spine for given msec
func (ss *Sim) NeuronUpdt(msec int, ge, gi float32, prespike bool) {
	ss.Msec = msec
	ly := ss.Net.LayerByName("Neuron").(axon.AxonLayer).AsAxon()
	nrn := ss.Neuron
	nex := &ss.NeuronEx

	vbio := chans.VToBio(nrn.Vm) // dend

	if prespike {
		ftime := float32(ss.Spine.States.Time)
		nex.PreSpike = 1
		if nex.PreSpikeT > 0 {
			nex.PreISI = ftime - nex.PreSpikeT
		}
		nex.PreSpikeT = ftime
	} else {
		nex.PreSpike = 0
	}

	// note: Ge should only
	nrn.GeRaw = ge
	ly.Act.Dt.GeSynFmRaw(nrn.GeRaw, &nrn.GeSyn, ly.Act.Init.Ge)
	nrn.Ge = nrn.GeSyn
	nrn.Gi = gi
	nrn.NMDA = ly.Act.NMDA.NMDA(nrn.NMDA, nrn.GeRaw, 1)
	nrn.Gnmda = ly.Act.NMDA.Gnmda(nrn.NMDA, nrn.VmDend) // gbar = 0 if !NMDAAxon
	nrn.GABAB, nrn.GABABx = ly.Act.GABAB.GABAB(nrn.GABAB, nrn.GABABx, nrn.Gi)
	nrn.GgabaB = ly.Act.GABAB.GgabaB(nrn.GABAB, nrn.VmDend)

	nex.Gvgcc = ss.VGCC.Gvgcc(nrn.VmDend, nex.VGCCm, nex.VGCCh)
	dm, dh := ss.VGCC.DMHFmV(nrn.VmDend, nex.VGCCm, nex.VGCCh)
	nex.VGCCm += dm
	nex.VGCCh += dh

	nex.Gak = ss.AK.Gak(nex.AKm, nex.AKh)
	dm, dh = ss.AK.DMHFmV(nrn.VmDend, nex.AKm, nex.AKh)
	nex.AKm += dm
	nex.AKh += dh

	nrn.Gk += nex.Gak
	nrn.Ge += nex.Gvgcc + nrn.Gnmda
	if !ss.NMDAAxon {
		nrn.Ge += ss.NMDAGbar * float32(ss.Spine.States.NMDAR.G)
	}
	nrn.Gi += nrn.GgabaB

	psd_pca := float32(1.7927e5 * 0.04) //  SVR_PSD
	cyt_pca := float32(1.0426e5 * 0.04) // SVR_CYT

	nex.VGCCJcaPSD = -vbio * psd_pca * nex.Gvgcc
	nex.VGCCJcaCyt = -vbio * cyt_pca * nex.Gvgcc

	ss.Spine.States.VmS = float64(vbio)

	ly.Act.VmFmG(nrn)
	ly.Act.ActFmG(nrn)

	// todo: Ca from NMDAAxon
	ss.Spine.Ca.SetInject(float64(nex.VGCCJcaPSD), float64(nex.VGCCJcaCyt))
	ss.Spine.States.PreSpike = float64(nex.PreSpike)

	ss.KinaseParams.Step(&ss.KinaseSyn, ss.Neuron, &ss.NeuronEx, float32(ss.Spine.States.CaSig.Ca.PSD))

	if !ss.KinaseOnly {
		ss.Spine.StepTime(0.001)
	}
}

// LogDefault does default logging for current Msec
func (ss *Sim) LogDefault(n int) {
	sfx := ""
	if n == 1 {
		sfx = "2"
	}
	msec := ss.Msec
	ss.LogTime(ss.Log("MsecLog"+sfx), msec%1000)
	if ss.Msec%10 == 0 {
		ss.LogTime(ss.Log("Msec10Log"+sfx), (msec/10)%1000)
		if ss.Msec%100 == 0 {
			ss.LogTime(ss.Log("Msec100Log"+sfx), (msec / 100))
			ss.UpdateTimePlots()
		}
	}
}

// Stop tells the sim to stop running
func (ss *Sim) Stop() {
	ss.StopNow = true
}

// Stopped is called when a run method stops running -- updates the IsRunning flag and toolbar
func (ss *Sim) Stopped() {
	ss.IsRunning = false
	if ss.Win != nil {
		vp := ss.Win.WinViewport2D()
		if ss.ToolBar != nil {
			ss.ToolBar.UpdateActions()
		}
		vp.SetNeedsFullRender()
	}
}

func (ss *Sim) GraphRun(secs float64, n int) {
	nms := int(secs / 0.001)
	sms := ss.Msec
	for msec := 0; msec < nms; msec++ {
		ss.NeuronUpdt(sms+msec, 0, 0, false)
		ss.LogDefault(n)
		if ss.StopNow {
			break
		}
	}
}

func (ss *Sim) RunQuiet(secs float64) {
	nms := int(secs / 0.001)
	sms := ss.Msec
	for msec := 0; msec < nms; msec++ {
		ss.NeuronUpdt(sms+msec, 0, 0, false)
		if ss.StopNow {
			break
		}
	}
}

//////////////////////////////////////////////
//  Time Log

// LogTime adds data from current msec to the given table at given row
func (ss *Sim) LogTime(dt *etable.Table, row int) {
	if dt.Rows <= row {
		dt.SetNumRows(row + 1)
	}
	nrn := ss.Neuron
	nex := &ss.NeuronEx

	dt.SetCellFloat("Time", row, float64(ss.Msec)*0.001)
	dt.SetCellFloat("Ge", row, float64(nrn.Ge))
	dt.SetCellFloat("Inet", row, float64(nrn.Inet))
	dt.SetCellFloat("Vm", row, float64(nrn.Vm))
	dt.SetCellFloat("Act", row, float64(nrn.Act))
	dt.SetCellFloat("Spike", row, float64(nrn.Spike))
	dt.SetCellFloat("Gk", row, float64(nrn.Gk))
	dt.SetCellFloat("ISI", row, float64(nrn.ISI))
	dt.SetCellFloat("AvgISI", row, float64(nrn.ISIAvg))
	dt.SetCellFloat("VmDend", row, float64(nrn.VmDend))
	dt.SetCellFloat("NMDA", row, float64(nrn.NMDA))
	dt.SetCellFloat("Gnmda", row, float64(nrn.Gnmda))
	dt.SetCellFloat("GABAB", row, float64(nrn.GABAB))
	dt.SetCellFloat("GgabaB", row, float64(nrn.GgabaB))
	dt.SetCellFloat("Gvgcc", row, float64(nex.Gvgcc))
	dt.SetCellFloat("VGCCm", row, float64(nex.VGCCm))
	dt.SetCellFloat("VGCCh", row, float64(nex.VGCCh))
	dt.SetCellFloat("VGCCJcaPSD", row, float64(nex.VGCCJcaPSD))
	dt.SetCellFloat("VGCCJcaCyt", row, float64(nex.VGCCJcaCyt))
	dt.SetCellFloat("Gak", row, float64(nex.Gak))
	dt.SetCellFloat("AKm", row, float64(nex.AKm))
	dt.SetCellFloat("AKh", row, float64(nex.AKh))

	ss.KinaseSyn.Log(dt, row)

	ss.Spine.Log(dt, row)
}

func (ss *Sim) ConfigTimeLog(dt *etable.Table) {
	dt.SetMetaData("name", "Kinase Time Log")
	dt.SetMetaData("desc", "Record of neuron / spine data over time")
	dt.SetMetaData("read-only", "true")
	dt.SetMetaData("precision", strconv.Itoa(LogPrec))

	sch := etable.Schema{
		{"Time", etensor.FLOAT64, nil, nil},
		{"Ge", etensor.FLOAT64, nil, nil},
		{"Inet", etensor.FLOAT64, nil, nil},
		{"Vm", etensor.FLOAT64, nil, nil},
		{"Act", etensor.FLOAT64, nil, nil},
		{"Spike", etensor.FLOAT64, nil, nil},
		{"Gk", etensor.FLOAT64, nil, nil},
		{"ISI", etensor.FLOAT64, nil, nil},
		{"AvgISI", etensor.FLOAT64, nil, nil},
		{"VmDend", etensor.FLOAT64, nil, nil},
		{"NMDA", etensor.FLOAT64, nil, nil},
		{"Gnmda", etensor.FLOAT64, nil, nil},
		{"GABAB", etensor.FLOAT64, nil, nil},
		{"GgabaB", etensor.FLOAT64, nil, nil},
		{"Gvgcc", etensor.FLOAT64, nil, nil},
		{"VGCCm", etensor.FLOAT64, nil, nil},
		{"VGCCh", etensor.FLOAT64, nil, nil},
		{"VGCCJcaPSD", etensor.FLOAT64, nil, nil},
		{"VGCCJcaCyt", etensor.FLOAT64, nil, nil},
		{"Gak", etensor.FLOAT64, nil, nil},
		{"AKm", etensor.FLOAT64, nil, nil},
		{"AKh", etensor.FLOAT64, nil, nil},
	}

	ss.KinaseSyn.ConfigLog(&sch)

	ss.Spine.ConfigLog(&sch)

	dt.SetFromSchema(sch, 1000)
}

func (ss *Sim) ConfigTimePlot(plt *eplot.Plot2D, dt *etable.Table) *eplot.Plot2D {
	plt.Params.Title = "Kinase Time Plot"
	plt.Params.XAxisCol = "Time"
	plt.SetTable(dt)
	// order of params: on, fixMin, min, fixMax, max
	plt.SetColParams("Time", eplot.Off, eplot.FloatMin, 0, eplot.FloatMax, 0)
	// plt.SetColParams("Time", eplot.Off, eplot.FixMin, .48, eplot.FixMax, .54)
	plt.SetColParams("Ge", eplot.Off, eplot.FixMin, 0, eplot.FixMax, 1)
	plt.SetColParams("Inet", eplot.Off, eplot.FixMin, -.2, eplot.FixMax, 1)
	plt.SetColParams("Vm", eplot.On, eplot.FixMin, 0, eplot.FixMax, 1)
	plt.SetColParams("Act", eplot.Off, eplot.FixMin, 0, eplot.FixMax, 1)
	plt.SetColParams("Spike", eplot.On, eplot.FixMin, 0, eplot.FixMax, 1)
	plt.SetColParams("Gk", eplot.Off, eplot.FixMin, 0, eplot.FixMax, 1)
	plt.SetColParams("ISI", eplot.Off, eplot.FixMin, -2, eplot.FloatMax, 1)
	plt.SetColParams("AvgISI", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("VmDend", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("NMDA", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("Gnmda", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("GABAB", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("GgabaB", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("Gvgcc", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("VGCCm", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("VGCCh", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("VGCCJcaPSD", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("VGCCJcaCyt", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("Gak", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("AKm", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("AKh", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)

	for _, cn := range dt.ColNames {
		if cn != "Time" {
			plt.SetColParams(cn, eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
		}
	}

	plt.SetColParams("VmS", eplot.Off, eplot.FixMin, -70, eplot.FloatMax, 1)
	plt.SetColParams("PreSpike", eplot.On, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("PSD_Ca", eplot.On, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("PSD_CaMact", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)

	plt.SetColParams("Ca", eplot.On, eplot.FixMin, 0, eplot.FloatMax, 2)
	plt.SetColParams("CaP", eplot.On, eplot.FixMin, 0, eplot.FloatMax, 2)
	plt.SetColParams("CaI", eplot.On, eplot.FixMin, 0, eplot.FloatMax, 2)
	plt.SetColParams("Ds", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 2)
	plt.SetColParams("Ps", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 2)
	plt.SetColParams("Wt", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 2)
	plt.SetColParams("DWt", eplot.Off, eplot.FloatMin, -1, eplot.FloatMax, 1)

	plt.SetColParams("Cyt_AC1act", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("PSD_AC1act", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("PSD_CaMKIIact", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("PSD_CaMKIIn2b", eplot.Off, eplot.FixMin, 0, eplot.FixMax, 10)
	if ss.Opts.UseDAPK1 {
		plt.SetColParams("PSD_DAPK1act", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
		plt.SetColParams("Cyt_DAPK1act", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
		plt.SetColParams("PSD_DAPK1n2b", eplot.Off, eplot.FixMin, 0, eplot.FixMax, 10)
	}
	plt.SetColParams("PSD_PP1act", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("PSD_CaNact", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("Cyt_CaMKIIact", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("Trp_AMPAR", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)

	return plt
}

func (ss *Sim) ResetTimePlots() {
	for nm, dt := range ss.Logs {
		if strings.HasPrefix(nm, "Msec") {
			dt.SetNumRows(0)
		}
	}
	ss.UpdateTimePlots()
}

func (ss *Sim) UpdateTimePlots() {
	for nm, plt := range ss.Plots {
		if strings.HasPrefix(nm, "Msec") {
			plt.Update()
		}
	}
}

//////////////////////////////////////////////
//  DWt Log

// LogDWt adds data for current dwt value as function of x, y values
func (ss *Sim) LogDWt(dt *etable.Table, x, y float64) {
	row := dt.Rows
	if dt.Rows <= row {
		dt.SetNumRows(row + 1)
	}

	dt.SetCellFloat("X", row, x)
	dt.SetCellFloat("Y", row, y)

	wt := ss.Spine.States.AMPAR.Trp.Tot
	dwt := (wt / ss.InitWt) - 1

	dt.SetCellFloat("DWt", row, float64(dwt))

	ss.Spine.Log(dt, row)
}

func (ss *Sim) ConfigDWtLog(dt *etable.Table) {
	dt.SetMetaData("name", "Kinase DWt Log")
	dt.SetMetaData("desc", "Record of final proportion dWt change")
	dt.SetMetaData("read-only", "true")
	dt.SetMetaData("precision", strconv.Itoa(LogPrec))

	sch := etable.Schema{
		{"X", etensor.FLOAT64, nil, nil},
		{"Y", etensor.FLOAT64, nil, nil},
		{"DWt", etensor.FLOAT64, nil, nil},
	}

	ss.Spine.ConfigLog(&sch)

	dt.SetFromSchema(sch, 0)
}

func (ss *Sim) ConfigDWtPlot(plt *eplot.Plot2D, dt *etable.Table) *eplot.Plot2D {
	plt.Params.Title = "Kinase DWt Plot"
	plt.Params.XAxisCol = "X"
	plt.Params.LegendCol = "Y"
	plt.SetTable(dt)
	// order of params: on, fixMin, min, fixMax, max
	plt.SetColParams("DWt", eplot.On, eplot.FixMin, -0.5, eplot.FixMax, 0.5)

	plt.SetColParams("PSD_CaMKIIact", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("PSD_PP1act", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("PSD_CaNact", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)

	return plt
}

//////////////////////////////////////////////
//  PhaseDWt Log

// LogPhaseDWt adds data for current dwt value as function of phase and phase hz levels
func (ss *Sim) LogPhaseDWt(dt *etable.Table, sphz, rphz []int) {
	row := dt.Rows
	if dt.Rows <= row {
		dt.SetNumRows(row + 1)
	}

	chl := (float64(sphz[1])/100.0)*(float64(rphz[1])/100.0) - (float64(sphz[0])/100.0)*(float64(rphz[0])/100.0)

	dt.SetCellFloat("CHL", row, float64(chl))
	dt.SetCellString("Cond", row, fmt.Sprintf("S:%d-%d R:%d-%d", sphz[0], sphz[1], rphz[0], rphz[1]))
	dt.SetCellFloat("SMhz", row, float64(sphz[0]))
	dt.SetCellFloat("SPhz", row, float64(sphz[1]))
	dt.SetCellFloat("RMhz", row, float64(rphz[0]))
	dt.SetCellFloat("RPhz", row, float64(rphz[1]))

	wt := ss.Spine.States.AMPAR.Trp.Tot
	dwt := (wt / ss.InitWt) - 1

	dt.SetCellFloat("DWt", row, float64(dwt))

	ss.Spine.Log(dt, row)
}

func (ss *Sim) ConfigPhaseDWtLog(dt *etable.Table) {
	dt.SetMetaData("name", "Kinase Phase DWt Log")
	dt.SetMetaData("desc", "Record of final proportion dWt change")
	dt.SetMetaData("read-only", "true")
	dt.SetMetaData("precision", strconv.Itoa(LogPrec))

	sch := etable.Schema{
		{"CHL", etensor.FLOAT64, nil, nil},
		{"Cond", etensor.STRING, nil, nil},
		{"SMhz", etensor.FLOAT64, nil, nil},
		{"SPhz", etensor.FLOAT64, nil, nil},
		{"RMhz", etensor.FLOAT64, nil, nil},
		{"RPhz", etensor.FLOAT64, nil, nil},
		{"DWt", etensor.FLOAT64, nil, nil},
	}

	ss.Spine.ConfigLog(&sch)

	dt.SetFromSchema(sch, 0)
}

func (ss *Sim) ConfigPhaseDWtPlot(plt *eplot.Plot2D, dt *etable.Table) *eplot.Plot2D {
	plt.Params.Title = "Kinase Phase DWt Plot"
	plt.Params.XAxisCol = "CHL"
	plt.Params.LegendCol = "Cond"
	plt.Params.Scale = 3
	plt.SetTable(dt)
	plt.Params.Points = true
	plt.Params.Lines = false
	// order of params: on, fixMin, min, fixMax, max
	plt.SetColParams("DWt", eplot.On, eplot.FixMin, -1, eplot.FixMax, 1)
	plt.SetColParams("CHL", eplot.Off, eplot.FixMin, -1, eplot.FixMax, 2)

	plt.SetColParams("PSD_CaMKIIact", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("PSD_PP1act", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)
	plt.SetColParams("PSD_CaNact", eplot.Off, eplot.FixMin, 0, eplot.FloatMax, 1)

	return plt
}

func (ss *Sim) ResetDWtPlot() {
	ss.Log("DWtLog").SetNumRows(0)
	ss.Plot("DWtPlot").Update()
	ss.Log("PhaseDWtLog").SetNumRows(0)
	ss.Plot("PhaseDWtPlot").Update()
}

//////////////////////////////////////////////
//  AutoK Log

func (ss *Sim) ConfigAutoKLog(dt *etable.Table) {
	dt.SetMetaData("name", "Kinase AutoK Plot")
	dt.SetMetaData("desc", "autoK as function of diff variables")
	dt.SetMetaData("read-only", "true")
	dt.SetMetaData("precision", strconv.Itoa(LogPrec))

	sch := etable.Schema{
		{"Total", etensor.FLOAT64, nil, nil},
		{"AutoK", etensor.FLOAT64, nil, nil},
	}

	dt.SetFromSchema(sch, 0)
}

func (ss *Sim) ConfigAutoKPlot(plt *eplot.Plot2D, dt *etable.Table) *eplot.Plot2D {
	plt.Params.Title = "Kinase AutoK Plot"
	plt.Params.XAxisCol = "Total"
	plt.SetTable(dt)
	// order of params: on, fixMin, min, fixMax, max
	plt.SetColParams("AutoK", eplot.On, eplot.FloatMin, 0, eplot.FloatMax, 1)
	return plt
}

// AutoK plots AutoK as function of T total
func (ss *Sim) AutoK() {
	dt := ss.Log("AutoKLog")

	max := 1.0
	inc := 0.01
	n := int(max / inc)
	dt.SetNumRows(n)

	cb := 0.75
	ct := 0.8
	ca := 0.8

	row := 0
	for T := 0.0; T <= max; T += inc {
		// this is a function that turns positive around 0.13 and is
		// roughly parabolic after that but not quite..
		tmp := T * (-0.22 + 1.826*T + -0.8*T*T)
		if false {
			Wb := .25 * T
			Wp := .25 * T
			Wt := .25 * T
			Wa := .25 * T
			tmp *= 0.75 * (cb*Wb + Wp + ct*Wt + ca*Wa)
			// if tmp < 0 {
			// 	tmp = 0
			// }
		}
		// autok := 0.29 * tmp
		autok := tmp
		dt.SetCellFloat("AutoK", row, autok)
		dt.SetCellFloat("Total", row, T)
		row++
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
// 		Gui

func (ss *Sim) ConfigNetView(nv *netview.NetView) {
	nv.ViewDefaults()
}

// ConfigGui configures the GoGi gui interface for this simulation,
func (ss *Sim) ConfigGui() *gi.Window {
	width := 1600
	height := 1200

	gi.SetAppName("kinase")
	gi.SetAppAbout(`This simulation replicates the Kinase et al, 2008 biophysical model of LTP / LTD.
See <a href="https://github.com/emer/axon/blob/master/examples/kinase/README.md">README.md on GitHub</a>.</p>`)

	win := gi.NewMainWindow("kinase", "Kinase", width, height)
	ss.Win = win

	vp := win.WinViewport2D()
	updt := vp.UpdateStart()

	mfr := win.SetMainFrame()

	tbar := gi.AddNewToolBar(mfr, "tbar")
	tbar.SetStretchMaxWidth()
	ss.ToolBar = tbar

	split := gi.AddNewSplitView(mfr, "split")
	split.Dim = mat32.X
	split.SetStretchMaxWidth()
	split.SetStretchMaxHeight()

	sv := giv.AddNewStructView(split, "sv")
	sv.SetStruct(ss)

	tv := gi.AddNewTabView(split, "tv")

	nv := tv.AddNewTab(netview.KiT_NetView, "NetView").(*netview.NetView)
	nv.Var = "Act"
	nv.SetNet(ss.Net)
	ss.NetView = nv
	ss.ConfigNetView(nv) // add labels etc

	plt := tv.AddNewTab(eplot.KiT_Plot2D, "DWtPlot").(*eplot.Plot2D)
	ss.AddPlot("DWtPlot", ss.ConfigDWtPlot(plt, ss.Log("DWtLog")))

	plt = tv.AddNewTab(eplot.KiT_Plot2D, "PhaseDWtPlot").(*eplot.Plot2D)
	ss.AddPlot("PhaseDWtPlot", ss.ConfigPhaseDWtPlot(plt, ss.Log("PhaseDWtLog")))

	plt = tv.AddNewTab(eplot.KiT_Plot2D, "Msec100Plot").(*eplot.Plot2D)
	ss.AddPlot("Msec100Plot", ss.ConfigTimePlot(plt, ss.Log("Msec100Log")))

	plt = tv.AddNewTab(eplot.KiT_Plot2D, "Msec100Plot2").(*eplot.Plot2D)
	ss.AddPlot("Msec100Plot2", ss.ConfigTimePlot(plt, ss.Log("Msec100Log2")))

	plt = tv.AddNewTab(eplot.KiT_Plot2D, "Msec10Plot").(*eplot.Plot2D)
	ss.AddPlot("Msec10Plot", ss.ConfigTimePlot(plt, ss.Log("Msec10Log")))

	plt = tv.AddNewTab(eplot.KiT_Plot2D, "Msec10Plot2").(*eplot.Plot2D)
	ss.AddPlot("Msec10Plot2", ss.ConfigTimePlot(plt, ss.Log("Msec10Log2")))

	plt = tv.AddNewTab(eplot.KiT_Plot2D, "MsecPlot").(*eplot.Plot2D)
	ss.AddPlot("MsecPlot", ss.ConfigTimePlot(plt, ss.Log("MsecLog")))

	plt = tv.AddNewTab(eplot.KiT_Plot2D, "MsecPlot2").(*eplot.Plot2D)
	ss.AddPlot("MsecPlot2", ss.ConfigTimePlot(plt, ss.Log("MsecLog2")))

	plt = tv.AddNewTab(eplot.KiT_Plot2D, "AutoKPlot").(*eplot.Plot2D)
	ss.AddPlot("AutoKPlot", ss.ConfigAutoKPlot(plt, ss.Log("AutoKLog")))

	split.SetSplits(.2, .8)

	tbar.AddAction(gi.ActOpts{Label: "Init", Icon: "update", Tooltip: "Initialize everything including network weights, and start over.  Also applies current params.", UpdateFunc: func(act *gi.Action) {
		act.SetActiveStateUpdt(!ss.IsRunning)
	}}, win.This(), func(recv, send ki.Ki, sig int64, data interface{}) {
		ss.Init()
		vp.SetNeedsFullRender()
	})

	tbar.AddAction(gi.ActOpts{Label: "Stop", Icon: "stop", Tooltip: "Interrupts running.  Hitting Train again will pick back up where it left off.", UpdateFunc: func(act *gi.Action) {
		act.SetActiveStateUpdt(ss.IsRunning)
	}}, win.This(), func(recv, send ki.Ki, sig int64, data interface{}) {
		ss.Stop()
	})

	tbar.AddAction(gi.ActOpts{Label: "Run", Icon: "step-fwd", Tooltip: "Runs current Stim.", UpdateFunc: func(act *gi.Action) {
		act.SetActiveStateUpdt(!ss.IsRunning)
	}}, win.This(), func(recv, send ki.Ki, sig int64, data interface{}) {
		if !ss.IsRunning {
			ss.IsRunning = true
			tbar.UpdateActions()
			ss.RunStim() // does go
		}
	})

	tbar.AddSeparator("run-sep")

	tbar.AddAction(gi.ActOpts{Label: "Reset Plots", Icon: "update", Tooltip: "Reset Time Plots.", UpdateFunc: func(act *gi.Action) {
		act.SetActiveStateUpdt(!ss.IsRunning)
	}}, win.This(), func(recv, send ki.Ki, sig int64, data interface{}) {
		if !ss.IsRunning {
			ss.ResetTimePlots()
		}
	})

	tbar.AddAction(gi.ActOpts{Label: "Reset DWt Plot", Icon: "update", Tooltip: "Reset DWt Plot.", UpdateFunc: func(act *gi.Action) {
		act.SetActiveStateUpdt(!ss.IsRunning)
	}}, win.This(), func(recv, send ki.Ki, sig int64, data interface{}) {
		if !ss.IsRunning {
			ss.ResetDWtPlot()
		}
	})

	tbar.AddAction(gi.ActOpts{Label: "AutoK Plot", Icon: "update", Tooltip: "Plot AutoK function.", UpdateFunc: func(act *gi.Action) {
		act.SetActiveStateUpdt(!ss.IsRunning)
	}}, win.This(), func(recv, send ki.Ki, sig int64, data interface{}) {
		if !ss.IsRunning {
			ss.AutoK()
		}
	})

	tbar.AddAction(gi.ActOpts{Label: "Defaults", Icon: "update", Tooltip: "Restore initial default parameters.", UpdateFunc: func(act *gi.Action) {
		act.SetActiveStateUpdt(!ss.IsRunning)
	}}, win.This(), func(recv, send ki.Ki, sig int64, data interface{}) {
		ss.Defaults()
		ss.Init()
		vp.SetNeedsFullRender()
	})

	tbar.AddAction(gi.ActOpts{Label: "README", Icon: "file-markdown", Tooltip: "Opens your browser on the README file that contains instructions for how to run this model."}, win.This(),
		func(recv, send ki.Ki, sig int64, data interface{}) {
			gi.OpenURL("https://github.com/emer/axon/blob/master/examples/kinase/README.md")
		})

	vp.UpdateEndNoSig(updt)

	// main menu
	appnm := gi.AppName()
	mmen := win.MainMenu
	mmen.ConfigMenus([]string{appnm, "File", "Edit", "Window"})

	amen := win.MainMenu.ChildByName(appnm, 0).(*gi.Action)
	amen.Menu.AddAppMenu(win)

	emen := win.MainMenu.ChildByName("Edit", 1).(*gi.Action)
	emen.Menu.AddCopyCutPaste(win)

	/*
		inQuitPrompt := false
		gi.SetQuitReqFunc(func() {
			if inQuitPrompt {
				return
			}
			inQuitPrompt = true
			gi.PromptDialog(vp, gi.DlgOpts{Title: "Really Quit?",
				Prompt: "Are you <i>sure</i> you want to quit and lose any unsaved params, weights, logs, etc?"}, gi.AddOk, gi.AddCancel,
				win.This(), func(recv, send ki.Ki, sig int64, data interface{}) {
					if sig == int64(gi.DialogAccepted) {
						gi.Quit()
					} else {
						inQuitPrompt = false
					}
				})
		})

		// gi.SetQuitCleanFunc(func() {
		// 	fmt.Printf("Doing final Quit cleanup here..\n")
		// })

		inClosePrompt := false
		win.SetCloseReqFunc(func(w *gi.Window) {
			if inClosePrompt {
				return
			}
			inClosePrompt = true
			gi.PromptDialog(vp, gi.DlgOpts{Title: "Really Close Window?",
				Prompt: "Are you <i>sure</i> you want to close the window?  This will Quit the App as well, losing all unsaved params, weights, logs, etc"}, gi.AddOk, gi.AddCancel,
				win.This(), func(recv, send ki.Ki, sig int64, data interface{}) {
					if sig == int64(gi.DialogAccepted) {
						gi.Quit()
					} else {
						inClosePrompt = false
					}
				})
		})
	*/

	win.SetCloseCleanFunc(func(w *gi.Window) {
		go gi.Quit() // once main window is closed, quit
	})

	win.MainMenuUpdated()
	return win
}

// These props register Save methods so they can be used
var SimProps = ki.Props{
	// "CallMethods": ki.PropSlice{
	// 	{"OpenGenesisData", ki.Props{
	// 		"desc": "Load data from Genesis version of Kinase model, from geneplot directory",
	// 		"icon": "file-open",
	// 		"Args": ki.PropSlice{
	// 			{"File Name", ki.Props{
	// 				"ext": ".tsv",
	// 			}},
	// 		},
	// 	}},
	// },
}

func mainrun() {
	TheSim.New()
	TheSim.Config()

	TheSim.Init()
	win := TheSim.ConfigGui()
	win.StartEventLoop()
}
