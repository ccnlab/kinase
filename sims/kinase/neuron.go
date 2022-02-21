// Copyright (c) 2021, The Emergent Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"

	"github.com/emer/axon/axon"
	"github.com/emer/axon/chans"
	"github.com/emer/emergent/chem"
	"github.com/emer/emergent/params"
)

// ParamSets for basic parameters
// Base is always applied, and others can be optionally selected to apply on top of that
var ParamSets = params.Sets{
	{Name: "Base", Desc: "these are the best params", Sheets: params.Sheets{
		"Network": &params.Sheet{
			{Sel: "Layer", Desc: "all defaults",
				Params: params.Params{
					"Layer.Act.Spike.Tr":     "7",
					"Layer.Act.Spike.RTau":   "3", // maybe could go a bit wider even
					"Layer.Act.NMDA.MgC":     "1.14",
					"Layer.Act.Decay.Glong":  "0.6", // 0.6
					"Layer.Act.Dend.GbarExp": "0.5", // 0.2 > 0.1 > 0
					"Layer.Act.Dend.GbarR":   "6",   // 3 > 2 good for 0.2 -- too low rel to ExpGbar causes fast ini learning, but then unravels
					"Layer.Act.Dt.VmDendTau": "5",   // 5 > 2.81 here but small effect
					"Layer.Act.Dt.GeTau":     "5",
					"Layer.Act.Dt.VmTau":     "1",
					"Layer.Act.Dt.VmSteps":   "2",
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
	NMDAGmg    float32 `desc:"NMDA mg-based blocking conductance"`
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
	LearnNow   float32 `desc:"when 0, it is time to learn according to theta cycle, otherwise increments up unless still -1 from init"`
}

func (nex *NeuronEx) Init() {
	nex.NMDAGmg = 0
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
	nex.LearnNow = -1
}

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
	nrn.GnmdaRaw = ge
	ly.Act.Dt.GeSynFmRaw(nrn.GeRaw, &nrn.GeSyn, ly.Act.Init.Ge)
	nrn.Ge = nrn.GeSyn
	nrn.Gi = gi
	ly.Act.NMDAFmRaw(nrn, 0)

	vmd := nrn.Vm
	if ss.DendVm {
		vmd = nrn.VmDend
	}

	nex.NMDAGmg = ly.Act.NMDA.MgGFmV(vmd)
	nrn.GABAB, nrn.GABABx = ly.Act.GABAB.GABAB(nrn.GABAB, nrn.GABABx, nrn.Gi)
	nrn.GgabaB = ly.Act.GABAB.GgabaB(nrn.GABAB, vmd)

	nex.Gvgcc = ss.VGCC.Gvgcc(vmd, nex.VGCCm, nex.VGCCh)
	dm, dh := ss.VGCC.DMHFmV(vmd, nex.VGCCm, nex.VGCCh)
	nex.VGCCm += dm
	nex.VGCCh += dh
	isi := nrn.ISI
	if isi >= ly.Act.Spike.VmR-1 && isi <= ly.Act.Spike.VmR {
		nex.VGCCm = 0 // resets
	}

	nex.Gak = ss.AK.Gak(nex.AKm, nex.AKh)
	dm, dh = ss.AK.DMHFmV(vmd, nex.AKm, nex.AKh)
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

	if !ss.KinaseOnly {
		ss.Spine.StepTime(0.001)
	}

	ss.KinaseParams.Step(&ss.KinaseSyn, ss.Neuron, &ss.NeuronEx, float32(chem.CoFmN(ss.Spine.States.CaSig.Ca.PSD, PSDVol)))
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
