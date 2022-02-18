// Copyright (c) 2021 The Emergent Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"github.com/emer/axon/axon"
	"github.com/emer/axon/chans"
	"github.com/emer/etable/etable"
	"github.com/emer/etable/etensor"
	"github.com/goki/ki/kit"
	"github.com/goki/mat32"
)

// KinaseState are synapse-specific Kinase algo state vars
type KinaseState struct {
	NMDAo float32 `desc:"Open NMDA from kinase"`
	NMDAi float32 `desc:"Inhibition of NMDA from presynaptic firing"`
	Ca    float32 `desc:"Computed Ca level"`
	Ds    float32 `desc:"LTD-driving DAPK1-based time integrator state"`
	Ps    float32 `desc:"LTP-driving CaMKII-based time integrator state"`
	CaM   float32 `desc:"Synaptic: first stage running average (mean) Ca calcium level (like CaM = calmodulin), feeds into CaP, for Kinase based learning"`
	CaP   float32 `desc:"Synaptic: shorter timescale integrated CaM value, representing the plus, LTP direction of weight change and capturing the function of CaMKII in the Kinase learning rule"`
	CaD   float32 `desc:"Synaptic: longer timescale integrated CaP value, representing the minus, LTD direction of weight change and capturing the function of DAPK1 in the Kinase learning rule"`
	RCaM  float32 `desc:"Recv: first stage running average (mean) Ca calcium level (like CaM = calmodulin), feeds into CaP, for Kinase based learning"`
	RCaP  float32 `desc:"Recv: shorter timescale integrated CaM value, representing the plus, LTP direction of weight change and capturing the function of CaMKII in the Kinase learning rule"`
	RCaD  float32 `desc:"Recv: longer timescale integrated CaP value, representing the minus, LTD direction of weight change and capturing the function of DAPK1 in the Kinase learning rule"`
	SCa   float32 `desc:"Computed Ca level - sending neuron"`
	SCaM  float32 `desc:"Sender: first stage running average (mean) Ca calcium level (like CaM = calmodulin), feeds into CaP, for Kinase based learning"`
	SCaP  float32 `desc:"Sender: shorter timescale integrated CaM value, representing the plus, LTP direction of weight change and capturing the function of CaMKII in the Kinase learning rule"`
	SCaD  float32 `desc:"Sender: longer timescale integrated CaP value, representing the minus, LTD direction of weight change and capturing the function of DAPK1 in the Kinase learning rule"`

	Wt  float32 `desc:"simulated weight"`
	DWt float32 `desc:"change in weight"`
}

func (ks *KinaseState) Init() {
	ks.Zero()
	ks.Wt = 0.5
}

func (ks *KinaseState) Zero() {
	ks.NMDAo = 0
	ks.NMDAi = 0
	ks.Ca = 0
	ks.Ds = 0
	ks.Ps = 0
	ks.CaM = 0
	ks.CaP = 0
	ks.CaD = 0
	ks.RCaM = 0
	ks.RCaP = 0
	ks.RCaD = 0
	ks.SCa = 0
	ks.SCaM = 0
	ks.SCaP = 0
	ks.SCaD = 0
	ks.Wt = 0
	ks.DWt = 0
}

func (ks *KinaseState) Log(dt *etable.Table, row int) {
	dt.SetCellFloat("NMDAo", row, float64(ks.NMDAo))
	dt.SetCellFloat("NMDAi", row, float64(ks.NMDAi))
	dt.SetCellFloat("Ca", row, float64(ks.Ca))
	dt.SetCellFloat("Ds", row, float64(ks.Ds))
	dt.SetCellFloat("Ps", row, float64(ks.Ps))
	dt.SetCellFloat("CaM", row, float64(ks.CaM))
	dt.SetCellFloat("CaP", row, float64(ks.CaP))
	dt.SetCellFloat("CaD", row, float64(ks.CaD))
	dt.SetCellFloat("RCaM", row, float64(ks.CaM))
	dt.SetCellFloat("RCaP", row, float64(ks.CaP))
	dt.SetCellFloat("RCaD", row, float64(ks.CaD))
	dt.SetCellFloat("SCa", row, float64(ks.Ca))
	dt.SetCellFloat("SCaM", row, float64(ks.CaM))
	dt.SetCellFloat("SCaP", row, float64(ks.CaP))
	dt.SetCellFloat("SCaD", row, float64(ks.CaD))
	if TheSim.KinaseLearn {
		dt.SetCellFloat("Wt", row, float64(ks.Wt))
		dt.SetCellFloat("DWt", row, float64(ks.DWt))
	}
}

func (ks *KinaseState) ConfigLog(sch *etable.Schema) {
	*sch = append(*sch, etable.Column{"NMDAo", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"NMDAi", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"Ca", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"Ds", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"Ps", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"CaM", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"CaP", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"CaD", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"RCaM", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"RCaP", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"RCaD", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"SCa", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"SCaM", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"SCaP", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"SCaD", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"Wt", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"DWt", etensor.FLOAT32, nil, nil})
}

/////////////////////////////////////////////////////////////////////

// KinaseRates are rate constants for integrating kinase states
type KinaseRates struct {
	Up float32 `desc:"rate of increase in state"`
	Dn float32 `desc:"rate of decrease / decay in state"`
}

func (kr *KinaseRates) Set(up, dn float32) {
	kr.Up = up
	kr.Dn = dn
}

// SetRR sets rate constants based on a rate factor and a ratio factor
// where ratio is how strong the up is relative to the down
func (kr *KinaseRates) SetRR(rate, ratio float32) {
	kr.Up = rate * ratio
	kr.Dn = rate
}

func (kr *KinaseRates) Step(s, cam float32) float32 {
	return 0.001 * (cam*kr.Up - s*kr.Dn)
}

// KinaseNMDA computes NMDA and resulting calcium from simulated allosteric NMDA receptors
// Targets the NMDAo Nopen number of open NMDA channels
type KinaseNMDA struct {
	CaGain   float32 `def:"2,5.5,300" desc:"overall Ca Gain factor: 2 for normalized with CaVdrive, 5.5 for CaVdrive, 300 without"`
	VGCCGain float32 `desc:"gain on Gvgcc for Ca influx"`
	PreOpen  float32 `def:"0.2" desc:"driver max in number open from presynaptic firing"`
	PreInhib float32 `def:"1" desc:"increment in inhibition from presynaptic firing"`
	DecayTau float32 `def:"30" desc:"conductance decay time constant"`
	InhibTau float32 `def:"100" desc:"presynaptic inhibition decay time constant"`
	CaVdrive bool    `desc:"when computing Ca, use the additional Vm-dependent drive factor -- reduces size of large spikes"`
}

func (kp *KinaseNMDA) Defaults() {
	kp.CaGain = 2 // 5.5 // 300
	kp.VGCCGain = 10
	kp.PreOpen = 0.2
	kp.PreInhib = 1
	kp.DecayTau = 30
	kp.InhibTau = 100
	kp.CaVdrive = true
}

// Vdrive returns the voltage-driven drive factor for computing Ca influx
func (kp *KinaseNMDA) Vdrive(vbio float32) float32 {
	if vbio > -0.1 && vbio < 0.1 {
		return (1.0 / (0.0756 + 0.5*vbio))
	} else {
		return -vbio / (1.0 - mat32.FastExp(0.0756*vbio))
	}
}

func (kp *KinaseNMDA) Step(ks *KinaseState, nrn *axon.Neuron, nex *NeuronEx) {
	vmd := nrn.Vm
	if TheSim.DendVm {
		vmd = nrn.VmDend
	}
	ks.NMDAo -= ks.NMDAo / kp.DecayTau
	ks.NMDAi -= ks.NMDAi / kp.InhibTau
	if nex.PreSpike > 0 {
		ks.NMDAo += (1 - ks.NMDAi) * (kp.PreOpen - ks.NMDAo)
		ks.NMDAi += (1 - ks.NMDAi) * kp.PreInhib
	}
	ks.Ca = kp.CaGain * ks.NMDAo * nex.NMDAGmg
	if kp.CaVdrive {
		ks.Ca *= kp.Vdrive(chans.VToBio(vmd))
	}
	ks.Ca += kp.VGCCGain * nex.Gvgcc
}

// KinaseRules are different options for Kinase-based learning rules
type KinaseRules int32

//go:generate stringer -type=KinaseRules

var KiT_KinaseRules = kit.Enums.AddEnum(KinaseRulesN, kit.NotBitFlag, nil)

func (ev KinaseRules) MarshalJSON() ([]byte, error)  { return kit.EnumMarshalJSON(ev) }
func (ev *KinaseRules) UnmarshalJSON(b []byte) error { return kit.EnumUnmarshalJSON(ev, b) }

// The time scales
const (
	// NeurSpkCa uses neuron-level spike-driven calcium signals
	// integrated at P vs. D time scales -- this is the original
	// Leabra and Axon XCAL / CHL learning rule.
	NeurSpkCa KinaseRules = iota

	// SynSpkCa uses synapse-level spike-driven calcium signals
	// with an OR rule for pre OR post spiking driving the CaM up,
	// which is then integrated at P vs. D time scales.
	// Basically a synapse version of original learning rule.
	SynSpkCa

	// SynSpkNMDA uses synapse-level spike-driven calcium signals
	// with an OR rule for pre OR post spiking driving the CaM up,
	// with NMDAo multiplying the spike drive to fit Bio Ca better
	// including the Bonus factor.
	// which is then integrated at P vs. D time scales.
	SynSpkNMDA

	// SynNMDACa uses synapse-level NMDA-driven calcium signals
	// (which can be either Urakubo allosteric or Kinase abstract)
	// integrated at P vs. D time scales -- abstract version
	// of the KinaseB biophysical learniung rule
	SynNMDACa

	KinaseRulesN
)

// KinaseSynParams has rate constants for averaging over activations
// at different time scales, to produce the running average activation
// values that then drive learning.
type KinaseSynParams struct {
	Rule    KinaseRules `desc:"which learning rule to use"`
	SpikeG  float32     `def:"20,8,200" desc:"spiking gain for Spk rules"`
	SBonus  float32     `def:"500" desc:"increment in spike impact as function of sender activity"`
	SAvgThr float32     `def:"0.02" desc:"optimization for compute speed -- threshold on sending avg values to update Ca values -- depends on Ca clearing upon Wt update"`
	MTau    float32     `def:"10,40" min:"1" desc:"CaM mean running-average time constant in cycles, which should be milliseconds typically (tau is roughly how long it takes for value to change significantly -- 1.4x the half-life). This provides a pre-integration step before integrating into the CaP short time scale"`
	PTau    float32     `def:"40,10" min:"1" desc:"LTP Ca-driven factor time constant in cycles, which should be milliseconds typically (tau is roughly how long it takes for value to change significantly -- 1.4x the half-life). Continuously updates based on current CaI value, resulting in faster tracking of plus-phase signals."`
	DTau    float32     `def:"40" min:"1" desc:"LTD Ca-driven factor time constant in cycles, which should be milliseconds typically (tau is roughly how long it takes for value to change significantly -- 1.4x the half-life).  Continuously updates based on current CaP value, resulting in slower integration that still reflects earlier minus-phase signals."`
	DScale  float32     `def:"1.05" desc:"scaling factor on CaD as it enters into the learning rule, to compensate for systematic decrease in activity over the course of a theta cycle"`

	MDt float32 `view:"-" json:"-" xml:"-" inactive:"+" desc:"rate = 1 / tau"`
	PDt float32 `view:"-" json:"-" xml:"-" inactive:"+" desc:"rate = 1 / tau"`
	DDt float32 `view:"-" json:"-" xml:"-" inactive:"+" desc:"rate = 1 / tau"`
}

func (kp *KinaseSynParams) Update() {
	kp.MDt = 1 / kp.MTau
	kp.PDt = 1 / kp.PTau
	kp.DDt = 1 / kp.DTau
}

func (kp *KinaseSynParams) Defaults() {
	kp.Rule = SynSpkCa
	kp.SpikeG = 8 // 200 // 8
	kp.SBonus = 500
	kp.SAvgThr = 0.02
	kp.MTau = 10 // 40
	kp.PTau = 40 // 10
	kp.DTau = 40
	kp.DScale = 1.05
	kp.Update()
}

// FmCa computes updates from current Ca value
func (kp *KinaseSynParams) FmCa(ca float32, caM, caP, caD *float32) {
	*caM += kp.MDt * (ca - *caM)
	*caP += kp.PDt * (*caM - *caP)
	*caD += kp.DDt * (*caP - *caD)
}

// DWt computes the weight change from CaP, CaD values
func (kp *KinaseSynParams) DWt(caP, caD float32) float32 {
	return caP - kp.DScale*caD
}

/////////////////////////////////////////////////////////////////////
// KinaseParams

// KinaseParams for abstract Kinase learning rule
type KinaseParams struct {
	NMDA  KinaseNMDA      `view:"inline" desc:"Ca computation params"`
	KinCa bool            `desc:"use the Kinase computed Ca instead of Urakubo"`
	SynCa KinaseSynParams `view:"inline" desc:"average Ca based on standard chained SS, S, M timescales"`
	Drate KinaseRates     `desc:"rates for Ds state for LTD / DAPK1 -- Ca versions"`
	Prate KinaseRates     `desc:"rates for Ps state for LTP / CaMKII -- Ca versions"`
	Lrate float32         `desc:"learning rate"`
}

func (kp *KinaseParams) Defaults() {
	kp.NMDA.Defaults()
	kp.KinCa = true
	kp.SynCa.Defaults()
	kp.Drate.Set(1.6, 0.7) // 1.6, 0.7
	kp.Prate.Set(1.8, 0.8) // matches CaMKII
	kp.Lrate = 0.02
}

// Step updates current state from params
func (kp *KinaseParams) Step(c *KinaseState, nrn *axon.Neuron, nex *NeuronEx, ca float32) {
	kp.NMDA.Step(c, nrn, nex)
	if kp.KinCa {
		ca = c.Ca
	}
	c.Ds += kp.Drate.Step(c.Ds, ca)
	c.Ps += kp.Prate.Step(c.Ps, ca)

	switch kp.SynCa.Rule {
	case NeurSpkCa:
		spk := kp.SynCa.SpikeG * nex.PreSpike
		kp.SynCa.FmCa(spk, &c.SCaM, &c.SCaP, &c.SCaD)
		spk = kp.SynCa.SpikeG * nrn.Spike
		kp.SynCa.FmCa(spk, &c.RCaM, &c.RCaP, &c.RCaD)
		c.CaM = c.RCaM * c.SCaM
		c.CaP = c.RCaP * c.SCaP
		c.CaD = c.RCaD * c.SCaD
	case SynSpkCa:
		spk := float32(0)
		if nrn.Spike > 0 || nex.PreSpike > 0 {
			spk = kp.SynCa.SpikeG
		}
		kp.SynCa.FmCa(spk, &c.CaM, &c.CaP, &c.CaD)
	case SynSpkNMDA:
		spk := float32(0)
		if nrn.Spike > 0 || nex.PreSpike > 0 {
			spk = kp.SynCa.SpikeG
			if nrn.Spike > 0 && nex.PreISI > 0 {
				spk += kp.SynCa.SBonus * c.NMDAo
			}
			spk *= c.NMDAo
		}
		kp.SynCa.FmCa(spk, &c.CaM, &c.CaP, &c.CaD)
	case SynNMDACa:
		c.CaM += kp.SynCa.MDt * (kp.SynCa.SpikeG*ca - c.CaM)
		c.CaP += kp.SynCa.PDt * (c.CaM - c.CaP)
		c.CaD += kp.SynCa.DDt * (c.CaP - c.CaD)
	}

	if nex.LearnNow == 0 {
		c.DWt = kp.Lrate * kp.SynCa.DWt(c.CaP, c.CaD) // todo: xcal
		c.Wt += c.DWt
	}
	if nex.LearnNow >= 0 {
		nex.LearnNow += 1
	}
}
