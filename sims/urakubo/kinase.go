// Copyright (c) 2021 The Emergent Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"github.com/emer/emergent/chem"
	"github.com/emer/etable/etable"
	"github.com/emer/etable/etensor"
)

type KinaseState struct {
	Ds float64 `desc:"LTD-driving DAPK1-based time integrator state"`
	Ps float64 `desc:"LTP-driving CaMKII-based time integrator state"`
	Wt float64 `desc:"simulated weight"`
}

func (ks *KinaseState) Init() {
	ks.Zero()
	ks.Wt = 1
}

func (ks *KinaseState) Zero() {
	ks.Ds = 0
	ks.Ps = 0
	ks.Wt = 0
}

func (ks *KinaseState) Integrate(d *KinaseState) {
	chem.Integrate(&ks.Ds, d.Ds)
	chem.Integrate(&ks.Ps, d.Ps)
}

func (ks *KinaseState) Log(dt *etable.Table, row int) {
	dt.SetCellFloat("Ds", row, ks.Ds)
	dt.SetCellFloat("Ps", row, ks.Ps)
	dt.SetCellFloat("Wt", row, ks.Wt)
}

func (ks *KinaseState) ConfigLog(sch *etable.Schema) {
	*sch = append(*sch, etable.Column{"Ds", etensor.FLOAT64, nil, nil})
	*sch = append(*sch, etable.Column{"Ps", etensor.FLOAT64, nil, nil})
	*sch = append(*sch, etable.Column{"Wt", etensor.FLOAT64, nil, nil})
}

// KinaseRates are rate constants for integrating kinase states
type KinaseRates struct {
	Up float64 `desc:"rate of increase in state"`
	Dn float64 `desc:"rate of decrease / decay in state"`
}

func (kr *KinaseRates) Set(up, dn float64) {
	kr.Up = up
	kr.Dn = dn
}

// SetRR sets rate constants based on a rate factor and a ratio factor
// where ratio is how strong the up is relative to the down
func (kr *KinaseRates) SetRR(rate, ratio float64) {
	kr.Up = rate * ratio
	kr.Dn = rate
}

func (kr *KinaseRates) Step(s, cam float64) float64 {
	return cam*kr.Up - s*kr.Dn
}

// KinaseParams for abstract Kinase learning rule
type KinaseParams struct {
	Drate KinaseRates `desc:"rates for Ds state for LTD / DAPK1"`
	Prate KinaseRates `desc:"rates for Ps state for LTP / CaMKII"`
	Lrate float64     `desc:"learning rate"`
}

func (kp *KinaseParams) Defaults() {
	kp.Drate.Set(1, 1) // CaM: 5, 5
	kp.Prate.Set(2, 3) // CaM: 8, 10
	kp.Lrate = 0.2
}

// Step computes deltas based on cam input
func (kp *KinaseParams) Step(c, d *KinaseState, cam, ca float64) {
	d.Ds = kp.Drate.Step(c.Ds, ca)
	d.Ps = kp.Prate.Step(c.Ps, ca)
}

// DWt computes weight change given values at this point
func (kp *KinaseParams) DWt(c *KinaseState) {
	c.Wt += kp.Lrate * (c.Ps - c.Ds)
}