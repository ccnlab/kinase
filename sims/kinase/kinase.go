// Copyright (c) 2021 The Emergent Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"github.com/emer/axon/axon"
	"github.com/emer/axon/chans"
	"github.com/emer/etable/etable"
	"github.com/emer/etable/etensor"
	"github.com/goki/mat32"
)

// KinaseState are synapse-specific Kinase algo state vars
type KinaseState struct {
	NMDAo   float32 `desc:"Open NMDA from kinase"`
	NMDAi   float32 `desc:"Inhibition of NMDA from presynaptic firing"`
	Ca      float32 `desc:"Computed Ca level"`
	Ds      float32 `desc:"LTD-driving DAPK1-based time integrator state"`
	Ps      float32 `desc:"LTP-driving CaMKII-based time integrator state"`
	AvgSS   float32 `desc:"super-short time-scale average of spiking -- goes up instantaneously when a Spike occurs, and then decays until the next spike -- provides the lowest-level time integration for running-averages that simulate accumulation of Calcium over time"`
	AvgS    float32 `desc:"short time-scale average of spiking, as a running average over AvgSS -- tracks the most recent activation states, and represents the plus phase for learning in error-driven learning (see AvgSLrn)"`
	AvgM    float32 `desc:"medium time-scale average of spiking, as a running average over AvgS -- represents the minus phase for error-driven learning"`
	AvgSLrn float32 `desc:"short time-scale activation average that is used for learning -- typically includes a small contribution from AvgMLrn in addition to mostly AvgS, as determined by LrnActAvgParams.LrnM -- important to ensure that when neuron turns off in plus phase (short time scale), enough medium-phase trace remains so that learning signal doesn't just go all the way to 0, at which point no learning would take place -- AvgS is subject to thresholding prior to mixing so low values become zero"`
	AvgMLrn float32 `desc:"medium time-scale activation average used in learning: subect to thresholding so low values become zero"`
	Wt      float32 `desc:"simulated weight"`
	DWt     float32 `desc:"change in weight"`
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
	ks.AvgSS = 0
	ks.AvgS = 0
	ks.AvgM = 0
	ks.AvgSLrn = 0
	ks.AvgMLrn = 0
	ks.Wt = 0
	ks.DWt = 0
}

func (ks *KinaseState) Log(dt *etable.Table, row int) {
	dt.SetCellFloat("NMDAo", row, float64(ks.NMDAo))
	dt.SetCellFloat("NMDAi", row, float64(ks.NMDAi))
	dt.SetCellFloat("Ca", row, float64(ks.Ca))
	dt.SetCellFloat("Ds", row, float64(ks.Ds))
	dt.SetCellFloat("Ps", row, float64(ks.Ps))
	dt.SetCellFloat("AvgSS", row, float64(ks.AvgSS))
	dt.SetCellFloat("AvgS", row, float64(ks.AvgS))
	dt.SetCellFloat("AvgM", row, float64(ks.AvgM))
	dt.SetCellFloat("AvgSLrn", row, float64(ks.AvgSLrn))
	dt.SetCellFloat("AvgMLrn", row, float64(ks.AvgMLrn))
	dt.SetCellFloat("Wt", row, float64(ks.Wt))
	dt.SetCellFloat("DWt", row, float64(ks.DWt))
}

func (ks *KinaseState) ConfigLog(sch *etable.Schema) {
	*sch = append(*sch, etable.Column{"NMDAo", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"NMDAi", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"Ca", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"Ds", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"Ps", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"AvgSS", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"AvgS", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"AvgM", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"AvgSLrn", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"AvgMLrn", etensor.FLOAT32, nil, nil})
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
	PreOpen  float32 `def:"0.25" desc:"driver max in number open from presynaptic firing"`
	PreInhib float32 `def:"1" desc:"increment in inhibition from presynaptic firing"`
	DecayTau float32 `def:"30" desc:"conductance decay time constant"`
	InhibTau float32 `def:"100" desc:"presynaptic inhibition decay time constant"`
	CaVdrive bool    `desc:"when computing Ca, use the additional Vm-dependent drive factor -- reduces size of large spikes"`
}

func (kp *KinaseNMDA) Defaults() {
	kp.CaGain = 2 // 5.5 // 300
	kp.PreOpen = 0.25
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
	ks.NMDAo -= ks.NMDAo / kp.DecayTau
	ks.NMDAi -= ks.NMDAi / kp.InhibTau
	if nex.PreSpike > 0 {
		ks.NMDAo += (1 - ks.NMDAi) * (kp.PreOpen - ks.NMDAo)
		ks.NMDAi += (1 - ks.NMDAi) * kp.PreInhib
	}
	ks.Ca = kp.CaGain * ks.NMDAo * nex.NMDAGmg
	if kp.CaVdrive {
		ks.Ca *= kp.Vdrive(chans.VToBio(nrn.VmDend))
	}
}

// LrnActAvgParams has rate constants for averaging over activations
// at different time scales, to produce the running average activation
// values that then drive learning in the XCAL learning rules.
// Is driven directly by spikes that increment running-average at super-short
// timescale.  Time cycle of 50 msec quarters / theta window learning works
// Cyc:50, SS:35 S:8, M:40 (best)
// Cyc:25, SS:20, S:4, M:20
type LrnActAvgParams struct {
	MinLrn float32 `def:"0.02" desc:"minimum learning activation -- below this goes to zero"`
	SSTau  float32 `def:"40" min:"1" desc:"time constant in cycles, which should be milliseconds typically (tau is roughly how long it takes for value to change significantly -- 1.4x the half-life), for continuously updating the super-short time-scale AvgSS value -- this is provides a pre-integration step before integrating into the AvgS short time scale -- it is particularly important for spiking -- in general 4 is the largest value without starting to impair learning, but a value of 7 can be combined with m_in_s = 0 with somewhat worse results"`
	STau   float32 `def:"10" min:"1" desc:"time constant in cycles, which should be milliseconds typically (tau is roughly how long it takes for value to change significantly -- 1.4x the half-life), for continuously updating the short time-scale AvgS value from the super-short AvgSS value (cascade mode) -- AvgS represents the plus phase learning signal that reflects the most recent past information"`
	MTau   float32 `def:"40" min:"1" desc:"time constant in cycles, which should be milliseconds typically (tau is roughly how long it takes for value to change significantly -- 1.4x the half-life), for continuously updating the medium time-scale AvgM value from the short AvgS value (cascade mode) -- AvgM represents the minus phase learning signal that reflects the expectation representation prior to experiencing the outcome (in addition to the outcome) -- the default value of 10 generally cannot be exceeded without impairing learning"`
	MScale float32 `min:"0" desc:"rescaling of M factor (multiplies AvgS when it drives M) to compensate for overall decrease in Ca influx over course of theta cycle"`
	LrnM   float32 `def:"0.1,0" min:"0" max:"1" desc:"how much of the medium term average activation to mix in with the short (plus phase) to compute the Neuron AvgSLrn variable that is used for the unit's short-term average in learning. This is important to ensure that when unit turns off in plus phase (short time scale), enough medium-phase trace remains so that learning signal doesn't just go all the way to 0, at which point no learning would take place -- typically need faster time constant for updating S such that this trace of the M signal is lost -- can set SSTau=7 and set this to 0 but learning is generally somewhat worse"`
	Init   float32 `def:"0.15" min:"0" max:"1" desc:"initial value for average"`

	SSDt float32 `view:"-" json:"-" xml:"-" inactive:"+" desc:"rate = 1 / tau"`
	SDt  float32 `view:"-" json:"-" xml:"-" inactive:"+" desc:"rate = 1 / tau"`
	MDt  float32 `view:"-" json:"-" xml:"-" inactive:"+" desc:"rate = 1 / tau"`
	LrnS float32 `view:"-" json:"-" xml:"-" inactive:"+" desc:"1-LrnM"`
}

// AvgsFmCa computes averages based on current Ca
func (aa *LrnActAvgParams) AvgsFmCa(ca float32, avgSS, avgS, avgM, avgSLrn, avgMLrn *float32) {
	*avgSS += aa.SSDt * (ca - *avgSS)
	*avgS += aa.SDt * (*avgSS - *avgS)
	*avgM += aa.MDt * (aa.MScale**avgS - *avgM)
	*avgMLrn = *avgM

	thrS := *avgS

	if *avgMLrn < aa.MinLrn && thrS < aa.MinLrn {
		*avgMLrn = 0
		thrS = 0
	}

	*avgSLrn = aa.LrnS*thrS + aa.LrnM**avgMLrn
}

func (aa *LrnActAvgParams) Update() {
	aa.SSDt = 1 / aa.SSTau
	aa.SDt = 1 / aa.STau
	aa.MDt = 1 / aa.MTau
	aa.LrnS = 1 - aa.LrnM
}

func (aa *LrnActAvgParams) Defaults() {
	aa.MinLrn = 0.02
	aa.SSTau = 40
	aa.STau = 10
	aa.MTau = 40
	aa.MScale = 0.93
	aa.LrnM = 0.1
	aa.Init = 0.15
	aa.Update()

}

/////////////////////////////////////////////////////////////////////
// KinaseParams

// KinaseParams for abstract Kinase learning rule
type KinaseParams struct {
	NMDA  KinaseNMDA      `view:"inline" desc:"Ca computation params"`
	KinCa bool            `desc:"use the Kinase computed Ca instead of Urakubo"`
	AvgCa LrnActAvgParams `view:"inline" desc:"average Ca based on standard chained SS, S, M timescales"`
	Drate KinaseRates     `desc:"rates for Ds state for LTD / DAPK1 -- Ca versions"`
	Prate KinaseRates     `desc:"rates for Ps state for LTP / CaMKII -- Ca versions"`
	Lrate float32         `desc:"learning rate"`
}

func (kp *KinaseParams) Defaults() {
	kp.NMDA.Defaults()
	kp.KinCa = true
	kp.AvgCa.Defaults()
	kp.Drate.Set(1.6, 0.7) // 1.6, 0.7
	kp.Prate.Set(1.8, 0.8) // matches CaMKII
	kp.Lrate = 0.2
}

// Step updates current state from params
func (kp *KinaseParams) Step(c *KinaseState, nrn *axon.Neuron, nex *NeuronEx, ca float32) {
	kp.NMDA.Step(c, nrn, nex)
	if kp.KinCa {
		ca = c.Ca
	}
	c.Ds += kp.Drate.Step(c.Ds, ca)
	c.Ps += kp.Prate.Step(c.Ps, ca)

	kp.AvgCa.AvgsFmCa(ca, &c.AvgSS, &c.AvgS, &c.AvgM, &c.AvgSLrn, &c.AvgMLrn)

	if nex.LearnNow == 0 {
		c.DWt = kp.Lrate * (c.AvgSLrn - c.AvgMLrn) // todo: xcal
		c.Wt += c.DWt
	}
	if nex.LearnNow >= 0 {
		nex.LearnNow += 1
	}
}
