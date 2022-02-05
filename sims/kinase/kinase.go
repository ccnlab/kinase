// Copyright (c) 2021 The Emergent Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"github.com/emer/axon/axon"
	"github.com/emer/etable/etable"
	"github.com/emer/etable/etensor"
)

// KinaseState are synapse-specific Kinase algo state vars
type KinaseState struct {
	Ca  float32 `desc:"Computed Ca level"`
	CaP float32 `desc:"Ca Pre binding"`
	CaI float32 `desc:"Ca Post inhibition"`
	Ds  float32 `desc:"LTD-driving DAPK1-based time integrator state"`
	Ps  float32 `desc:"LTP-driving CaMKII-based time integrator state"`
	Wt  float32 `desc:"simulated weight"`
	DWt float32 `desc:"change in weight"`
}

func (ks *KinaseState) Init() {
	ks.Zero()
	ks.Wt = 1
}

func (ks *KinaseState) Zero() {
	ks.Ca = 0
	ks.CaP = 0
	ks.CaI = 0
	ks.Ds = 0
	ks.Ps = 0
	ks.Wt = 0
	ks.DWt = 0
}

func (ks *KinaseState) Log(dt *etable.Table, row int) {
	dt.SetCellFloat("Ca", row, float64(ks.Ca))
	dt.SetCellFloat("CaP", row, float64(ks.CaP))
	dt.SetCellFloat("CaI", row, float64(ks.CaI))
	dt.SetCellFloat("Ds", row, float64(ks.Ds))
	dt.SetCellFloat("Ps", row, float64(ks.Ps))
	dt.SetCellFloat("Wt", row, float64(ks.Wt))
	dt.SetCellFloat("DWt", row, float64(ks.DWt))
}

func (ks *KinaseState) ConfigLog(sch *etable.Schema) {
	*sch = append(*sch, etable.Column{"Ca", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"CaP", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"CaI", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"Ds", etensor.FLOAT32, nil, nil})
	*sch = append(*sch, etable.Column{"Ps", etensor.FLOAT32, nil, nil})
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

// KinaseCa computes calcium from simulated allosteric NMDA receptors
type KinaseCa struct {
	PreUp          float32 `desc:"increment from presynaptic firing"`
	PostUp         float32 `desc:"postsynaptic firing increment by itself"`
	PostPreUp      float32 `desc:"extra bonus from residual pre"`
	InhibTau       float32 `desc:"postsynaptic inhibition time constant"`
	PostInhibDelay float32 `desc:"delay in msec for onset of post inhibition"`
	DecayTau       float32 `desc:"decay time constant"`
	PreDecayTau    float32 `desc:"rate of presynaptic decay"`
}

func (kc *KinaseCa) Defaults() {
	kc.PreUp = 4.5
	kc.PostUp = 0.9
	kc.PostPreUp = 15
	kc.InhibTau = 50
	kc.PostInhibDelay = 4
	kc.DecayTau = 5
	kc.PreDecayTau = 30
}

func (kc *KinaseCa) Step(ks *KinaseState, nrn *axon.Neuron, nex *NeuronEx) {
	ifact := (1 - ks.CaI)
	dkt := kc.DecayTau + ks.CaP*kc.PreDecayTau

	ks.Ca -= (ks.Ca - ks.CaP) / dkt
	ks.CaI -= ks.CaI / kc.InhibTau
	ks.CaP -= ks.CaP / kc.PreDecayTau
	if nex.PreSpike > 0 {
		if ks.CaP > 0.1 {
			ks.CaI = 1
		}
		ks.Ca += ifact * (1 - ks.CaP) * kc.PreUp
		ks.CaP = 1
	}
	if nrn.Spike > 0 { // post
		ks.Ca += kc.PostUp + ifact*kc.PostPreUp*ks.CaP
		ks.CaP = 0 // resets
	}
	if nrn.ISI == kc.PostInhibDelay {
		ks.CaI = 1
	}
}

/////////////////////////////////////////////////////////////////////
// KinaseParams

// KinaseParams for abstract Kinase learning rule
type KinaseParams struct {
	Ca    KinaseCa    `view:"inline" desc:"Ca computation params"`
	Drate KinaseRates `desc:"rates for Ds state for LTD / DAPK1 -- Ca versions"`
	Prate KinaseRates `desc:"rates for Ps state for LTP / CaMKII -- Ca versions"`
	Lrate float32     `desc:"learning rate"`
}

func (kp *KinaseParams) Defaults() {
	kp.Ca.Defaults()
	kp.Drate.Set(0.55, 0.35) // Ca matching CaMKII = 0.55, 0.35
	kp.Prate.Set(0.6, 0.4)   // Ca matching CaMKII = 0.6, 0.4
	kp.Lrate = 0.2
}

// Step updates current state from params
func (kp *KinaseParams) Step(c *KinaseState, nrn *axon.Neuron, nex *NeuronEx, ca float32) {
	kp.Ca.Step(c, nrn, nex)
	c.Ds += kp.Drate.Step(c.Ds, ca)
	c.Ps += kp.Prate.Step(c.Ps, ca)

	// todo:
	// if pi == 1 && (dur-msec) == 50 {
	// 	ss.Spine.Kinase.DWt(&ss.Spine.States.Kinase)
	// }
}

// DWt computes weight change given values at this point
func (kp *KinaseParams) DWt(c *KinaseState) {
	c.DWt = kp.Lrate * (c.Ps - c.Ds)
	c.Wt += c.DWt
}
