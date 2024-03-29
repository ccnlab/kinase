// Copyright (c) 2021, The Emergent Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import "github.com/emer/emergent/params"

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
					"Layer.Act.AK.Gbar":   "0",
				}},
		},
	}},
}

// Extra state for neuron -- Vgcc and AK
type NeuronEx struct {
	VgccJcaPSD float32 `desc:"Vgcc Ca calcium contribution to PSD"`
	VgccJcaCyt float32 `desc:"Vgcc Ca calcium contribution to Cyt"`
	AKm        float32 `desc:"AK M gate -- activates with increasing Vm"`
	AKh        float32 `desc:"AK H gate -- deactivates with increasing Vm"`
}

func (nex *NeuronEx) Init() {
	nex.VgccJcaPSD = 0
	nex.VgccJcaCyt = 0
	nex.AKm = 0
	nex.AKh = 1
}
