// Copyright (c) 2021 The Emergent Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"math/rand"

	"github.com/goki/ki/kit"
	"github.com/goki/mat32"
)

type Stims int32

//go:generate stringer -type=Stims

var KiT_Stims = kit.Enums.AddEnum(StimsN, kit.NotBitFlag, nil)

func (ev Stims) MarshalJSON() ([]byte, error)  { return kit.EnumMarshalJSON(ev) }
func (ev *Stims) UnmarshalJSON(b []byte) error { return kit.EnumUnmarshalJSON(ev, b) }

// The different stimulus functions
const (
	Baseline Stims = iota

	STDP

	STDPSweep

	STDPPacketSweep

	Poisson

	SPoissonRGClamp

	PoissonHzSweep

	PoissonDurSweep

	OpPhaseDurSweep

	ThetaErr

	ThetaErrComp

	ThetaErrSweep

	ThetaErrAllSweep

	StimsN
)

// StimFuncs are the stimulus functions
var StimFuncs = map[Stims]func(){
	Baseline:         BaselineFun,
	STDP:             STDPFun,
	STDPSweep:        STDPSweepFun,
	STDPPacketSweep:  STDPPacketSweepFun,
	Poisson:          PoissonFun,
	SPoissonRGClamp:  SPoissonRGClampFun,
	PoissonHzSweep:   PoissonHzSweepFun,
	PoissonDurSweep:  PoissonDurSweepFun,
	OpPhaseDurSweep:  OpPhaseDurSweepFun,
	ThetaErr:         ThetaErrFun,
	ThetaErrComp:     ThetaErrCompFun,
	ThetaErrSweep:    ThetaErrSweepFun,
	ThetaErrAllSweep: ThetaErrAllSweepFun,
}

// RGeStimForHzMap is the strength of GeStim G clamp to obtain a given R firing rate
var RGeStimForHzMap = map[int]float32{
	25:  .09,
	50:  .12,
	100: .15,
}

func RGeStimForHz(hz float32) float32 {
	var gel, geh, hzl, hzh float32
	switch {
	case hz <= 25:
		gel = 0
		geh = RGeStimForHzMap[25]
		hzl = 0
		hzh = 25
	case hz <= 50:
		gel = RGeStimForHzMap[25]
		geh = RGeStimForHzMap[50]
		hzl = 25
		hzh = 50
	case hz <= 100:
		gel = RGeStimForHzMap[50]
		geh = RGeStimForHzMap[100]
		hzl = 50
		hzh = 100
	default:
		gel = RGeStimForHzMap[100]
		geh = 2 * gel
		hzl = 100
		hzh = 200
	}
	return gel + ((hz-hzl)/(hzh-hzl))*(geh-gel)
}

func BaselineFun() {
	ss := &TheSim
	for msec := 0; msec < 500000; msec++ { // 500000 = 500 sec for full baseline
		ss.NeuronUpdt(msec, 0, 0, false)
		ss.LogDefault(0)
		if ss.StopNow {
			break
		}
	}
	ss.Spine.InitCode()
	ss.Stopped()
}

func STDPFun() {
	ss := &TheSim
	toff := 500
	dur := 1
	psms := toff + 5 - ss.DeltaT // 5 is lag
	tott := ss.NReps * 1000

	for msec := 0; msec < tott; msec++ {
		ims := msec % 1000
		prespike := ims == psms
		ge := float32(0.0)
		if ims >= toff && ims < toff+dur {
			ge = ss.GeStim
		}
		if ims == toff+ss.STDPLrnOff {
			ss.LearnNow()
		}
		ss.NeuronUpdt(msec, ge, 0, prespike)
		ss.LogDefault(0)
		if ss.StopNow {
			break
		}
	}
	ss.GraphRun(ss.FinalSecs, 0)
	ss.Stopped()
}

func STDPSweepFun() {
	ss := &TheSim
	toff := 500
	dur := 1
	tott := ss.NReps * 1000

	ss.ResetDWtPlot()

	for dt := -ss.DeltaTRange; dt <= ss.DeltaTRange; dt += ss.DeltaTInc {
		psms := toff + 5 - dt // 5 is lag
		ss.ResetTimePlots()
		ss.Init()

		for msec := 0; msec < tott; msec++ {
			ims := msec % 1000
			prespike := ims == psms
			ge := float32(0.0)
			if ims >= toff && ims < toff+dur {
				ge = ss.GeStim
			}
			if ims == toff+ss.STDPLrnOff {
				ss.LearnNow()
			}
			ss.NeuronUpdt(msec, ge, 0, prespike)
			ss.LogDefault(0)
			if ss.StopNow {
				ss.Stopped()
				return
			}
		}
		ss.GraphRun(ss.FinalSecs, 0)
		ss.LogDWt(ss.Log("DWtLog"), float64(dt), 0)
		ss.Plot("DWtPlot").GoUpdate()
	}

	ss.Stopped()
}

// STDPPacket runs a sequence of Dur pre-post spike packets with sweep of
// pre-post offset in < 1/2 SendHz ISI range, with ISI interval between packets, N reps,
// and varying the frequency of pre-post firing (X axis).
func STDPPacketSweepFun() {
	ss := &TheSim

	isi := int(1000.0 / ss.SendHz)
	hisi := isi / 2
	dr := hisi - 5 // allow for lag

	ss.ResetDWtPlot()

	for dt := -dr; dt <= dr; dt++ {
		rms := hisi
		sms := hisi + 5 - dt // 5 is lag
		ss.ResetTimePlots()
		ss.Init()

		for ri := 0; ri < ss.NReps; ri++ {
			for msec := 0; msec < ss.DurMsec; msec++ {
				ims := msec % isi
				prespike := ims == sms
				ge := float32(0.0)
				if ims == rms {
					ge = ss.GeStim
				}
				ss.NeuronUpdt(msec, ge, 0, prespike)
				ss.LogDefault(0)
				if ss.StopNow {
					ss.Stopped()
					return
				}
			}
			ss.LearnNow()
		}
		ss.GraphRun(ss.FinalSecs, 0)
		ss.LogDWt(ss.Log("DWtLog"), float64(dt), float64(ss.SendHz))
		ss.Plot("DWtPlot").GoUpdate()
	}

	ss.Stopped()
}

func PoissonFun() {
	ss := &TheSim

	Sint := mat32.Exp(-1000.0 / ss.SendHz)
	Rint := mat32.Exp(-1000.0 / ss.RecvHz)

	tmsec := 0
	for ri := 0; ri < ss.NReps; ri++ {
		Sp := float32(1)
		Rp := float32(1)

		for msec := 0; msec < ss.DurMsec; msec++ {
			Sp *= rand.Float32()
			prespike := false
			if Sp <= Sint {
				prespike = true
				Sp = 1
			}

			ge := float32(0.0)
			Rp *= rand.Float32()
			if Rp <= Rint {
				ge = ss.GeStim
				Rp = 1
			}

			ss.NeuronUpdt(tmsec, ge, 0, prespike)
			ss.LogDefault(0)
			if ss.StopNow {
				break
			}
			tmsec++
		}
		ss.LearnNow()
		ss.GraphRun(ss.ISISec, 0)
	}
	ss.GraphRun(ss.FinalSecs, 0)
	ss.Stopped()
}

func SPoissonRGClampFun() {
	ss := &TheSim

	Sint := mat32.Exp(-1000.0 / ss.SendHz)

	for ri := 0; ri < ss.NReps; ri++ {
		Sp := float32(1)

		for msec := 0; msec < ss.DurMsec; msec++ {
			Sp *= rand.Float32()
			prespike := false
			if Sp <= Sint {
				prespike = true
				Sp = 1
			}

			ss.NeuronUpdt(msec, ss.GeStim, 0, prespike)
			ss.LogDefault(0)
			if ss.StopNow {
				break
			}
		}
		ss.LearnNow()
		ss.GraphRun(ss.ISISec, 0)
	}
	ss.GraphRun(ss.FinalSecs, 0)
	ss.Stopped()
}

func PoissonHzSweepFun() {
	ss := &TheSim

	ss.ResetDWtPlot()

	for shz := 10; shz <= 100; shz += 10 {
		for rhz := 10; rhz <= 100; rhz += 10 {
			Sint := mat32.Exp(-1000.0 / float32(shz))
			Rint := mat32.Exp(-1000.0 / float32(rhz))

			ss.ResetTimePlots()
			ss.Init()
			for ri := 0; ri < ss.NReps; ri++ {
				Sp := float32(1)
				Rp := float32(1)

				for msec := 0; msec < ss.DurMsec; msec++ {
					Sp *= rand.Float32()
					prespike := false
					if Sp <= Sint {
						prespike = true
						Sp = 1
					}

					ge := float32(0.0)
					Rp *= rand.Float32()
					if Rp <= Rint {
						ge = ss.GeStim
						Rp = 1
					}

					ss.NeuronUpdt(msec, ge, 0, prespike)
					ss.LogDefault(0)
					if ss.StopNow {
						ss.Stopped()
						return
					}
				}
				ss.LearnNow()
				ss.GraphRun(ss.ISISec, 0)
			}
			ss.GraphRun(ss.FinalSecs, 0)
			ss.LogDWt(ss.Log("DWtLog"), float64(rhz), float64(shz))
			ss.Plot("DWtPlot").GoUpdate()
		}
	}
	ss.Stopped()
}

func PoissonDurSweepFun() {
	ss := &TheSim

	ss.ResetDWtPlot()

	for dur := 200; dur <= 1000; dur += 100 {
		for rhz := 10; rhz <= 100; rhz += 10 {
			Sint := mat32.Exp(-1000.0 / float32(ss.SendHz))
			Rint := mat32.Exp(-1000.0 / float32(rhz))

			ss.ResetTimePlots()
			ss.Init()
			for ri := 0; ri < ss.NReps; ri++ {
				Sp := float32(1)
				Rp := float32(1)

				for msec := 0; msec < dur; msec++ {
					Sp *= rand.Float32()
					prespike := false
					if Sp <= Sint {
						prespike = true
						Sp = 1
					}

					ge := float32(0.0)
					Rp *= rand.Float32()
					if Rp <= Rint {
						ge = ss.GeStim
						Rp = 1
					}

					ss.NeuronUpdt(msec, ge, 0, prespike)
					ss.LogDefault(0)
					if ss.StopNow {
						ss.Stopped()
						return
					}
				}
				ss.LearnNow()
				ss.GraphRun(ss.ISISec, 0)
			}
			ss.GraphRun(ss.FinalSecs, 0)
			ss.LogDWt(ss.Log("DWtLog"), float64(rhz), float64(dur))
			ss.Plot("DWtPlot").GoUpdate()
		}
	}
	ss.Stopped()
}

// OpPhase runs sending, recv in opposite phases (half interval off at start)
// This is what was used in the original XCAL Dwt function derivation in Genesis model
func OpPhaseDurSweepFun() {
	ss := &TheSim

	ss.ResetDWtPlot()

	for dur := 200; dur <= 1000; dur += 100 {
		for rhz := 10; rhz <= 100; rhz += 10 {
			Sint := 1000.0 / float32(ss.SendHz)
			Rint := 1000.0 / float32(rhz)

			ss.ResetTimePlots()
			ss.Init()
			for ri := 0; ri < ss.NReps; ri++ {
				Sp := Sint / 2
				Rp := Rint

				for msec := 0; msec < dur; msec++ {
					fms := float32(msec)
					prespike := false
					if fms-Sp >= Sint {
						prespike = true
						Sp = fms
					}

					ge := float32(0.0)
					if fms-Rp >= Rint {
						ge = ss.GeStim
						Rp = fms
					}

					ss.NeuronUpdt(msec, ge, 0, prespike)
					ss.LogDefault(0)
					if ss.StopNow {
						ss.Stopped()
						return
					}
				}
				ss.LearnNow()
				ss.GraphRun(ss.ISISec, 0)
			}
			ss.GraphRun(ss.FinalSecs, 0)
			ss.LogDWt(ss.Log("DWtLog"), float64(rhz), float64(dur))
			ss.Plot("DWtPlot").GoUpdate()
		}
	}
	ss.Stopped()
}

func ThetaErrFun() {
	ss := &TheSim

	phsdur := []int{ss.DurMsec / 2, ss.DurMsec / 2}
	nphs := len(phsdur)

	// using send, recv for minus, plus
	sphz := []int{int(ss.SendHz), int(ss.RecvHz)}
	rphz := []int{int(ss.SendHz), int(ss.RecvHz)}

	tmsec := 0
	ss.ResetTimePlots()
	ss.Init()
	ss.RunQuiet(10)
	for ri := 0; ri < ss.NReps; ri++ {
		Sp := float32(1)
		Rp := float32(1)
		for pi := 0; pi < nphs; pi++ {
			dur := phsdur[pi]
			shz := sphz[pi]
			rhz := rphz[pi]
			Sint := mat32.Exp(-1000.0 / float32(shz))
			Rint := mat32.Exp(-1000.0 / float32(rhz))
			for msec := 0; msec < dur; msec++ {
				Sp *= rand.Float32()
				prespike := false
				if Sp <= Sint {
					prespike = true
					Sp = 1
				}

				ge := float32(0.0)
				Rp *= rand.Float32()
				if Rp <= Rint {
					ge = ss.GeStim
					Rp = 1
				}
				if ss.RGClamp {
					ge = RGeStimForHz(float32(rhz))
				}

				ss.NeuronUpdt(tmsec, ge, 0, prespike)
				ss.LogDefault(0)
				if ss.StopNow {
					ss.Stopped()
					return
				}
				tmsec++
			}
		}
		ss.LearnNow()
		ss.GraphRun(ss.ISISec, 0)
		tmsec = ss.Msec
	}
	ss.GraphRun(ss.FinalSecs, 0)
	tmsec = ss.Msec
	ss.Stopped()
}

func ThetaErrCompFun() {
	ss := &TheSim

	ss.ResetTimePlots()
	for itr := 0; itr < 2; itr++ {
		phsdur := []int{ss.DurMsec / 2, ss.DurMsec / 2}
		nphs := len(phsdur)

		// using send, recv for minus, plus
		sphz := []int{int(ss.SendHz), int(ss.RecvHz)}
		rphz := []int{int(ss.SendHz), int(ss.RecvHz)}

		if itr == 1 {
			sphz[1] = int(ss.SendHz)
			rphz[1] = int(ss.SendHz)
		}

		tmsec := 0
		ss.Init()
		// ss.RunQuiet(10)

		for ri := 0; ri < ss.NReps; ri++ {
			Sp := float32(1)
			Rp := float32(1)
			for pi := 0; pi < nphs; pi++ {
				dur := phsdur[pi]
				shz := sphz[pi]
				rhz := rphz[pi]
				Sint := mat32.Exp(-1000.0 / float32(shz))
				Rint := mat32.Exp(-1000.0 / float32(rhz))
				for msec := 0; msec < dur; msec++ {
					Sp *= rand.Float32()
					prespike := false
					if Sp <= Sint {
						prespike = true
						Sp = 1
					}

					ge := float32(0.0)
					Rp *= rand.Float32()
					if Rp <= Rint {
						ge = ss.GeStim
						Rp = 1
					}
					if ss.RGClamp {
						ge = RGeStimForHz(float32(rhz))
					}

					ss.NeuronUpdt(tmsec, ge, 0, prespike)
					ss.LogDefault(itr)
					if ss.StopNow {
						ss.Stopped()
						return
					}
					tmsec++
				}
			}
			ss.LearnNow()
			ss.GraphRun(ss.ISISec, itr)
			tmsec = ss.Msec
		}
		ss.GraphRun(ss.FinalSecs, itr)
		tmsec = ss.Msec
		ss.LogDWt(ss.Log("DWtLog"), float64(itr), 0)
		ss.Plot("DWtPlot").GoUpdate()
	}
	ss.Stopped()
}

func ThetaErrSweepFun() {
	ss := &TheSim

	ss.ResetDWtPlot()

	hz := []int{25, 50, 100}
	nhz := len(hz)

	phsdur := []int{ss.DurMsec / 2, ss.DurMsec / 2}
	nphs := len(phsdur)

	sphz := []int{0, 0}
	rphz := []int{0, 0}

	for smi := 0; smi < nhz; smi++ {
		sphz[0] = hz[smi] // minus phase
		rphz[0] = hz[smi] // minus phase
		for spi := 0; spi < nhz; spi++ {
			sphz[1] = hz[spi] // plus phase
			rphz[1] = hz[spi] // plus phase

			tmsec := 0
			ss.ResetTimePlots()
			ss.Init()
			for ri := 0; ri < ss.NReps; ri++ {
				Sp := float32(1)
				Rp := float32(1)
				for pi := 0; pi < nphs; pi++ {
					dur := phsdur[pi]
					shz := sphz[pi]
					rhz := rphz[pi]
					Sint := mat32.Exp(-1000.0 / float32(shz))
					Rint := mat32.Exp(-1000.0 / float32(rhz))
					for msec := 0; msec < dur; msec++ {
						Sp *= rand.Float32()
						prespike := false
						if Sp <= Sint {
							prespike = true
							Sp = 1
						}

						ge := float32(0.0)
						Rp *= rand.Float32()
						if Rp <= Rint {
							ge = ss.GeStim
							Rp = 1
						}
						if ss.RGClamp {
							ge = RGeStimForHz(float32(rhz))
						}

						ss.NeuronUpdt(tmsec, ge, 0, prespike)
						ss.LogDefault(0)
						if ss.StopNow {
							ss.Stopped()
							return
						}
						tmsec++
					}
				}
				ss.LearnNow()
				ss.GraphRun(ss.ISISec, 0)
				tmsec = ss.Msec
			}
			ss.GraphRun(ss.FinalSecs, 0)
			tmsec = ss.Msec
			ss.LogPhaseDWt(ss.Log("PhaseDWtLog"), sphz, rphz)
			ss.Plot("PhaseDWtPlot").GoUpdate()
		}
	}
	ss.Stopped()
}

func ThetaErrAllSweepFun() {
	ss := &TheSim

	ss.ResetDWtPlot()

	hz := []int{25, 50, 100}
	nhz := len(hz)

	phsdur := []int{ss.DurMsec / 2, ss.DurMsec / 2}
	nphs := len(phsdur)

	sphz := []int{0, 0}
	rphz := []int{0, 0}

	for smi := 0; smi < nhz; smi++ {
		sphz[0] = hz[smi] // minus phase
		for spi := 0; spi < nhz; spi++ {
			sphz[1] = hz[spi] // plus phase
			for rmi := 0; rmi < nhz; rmi++ {
				rphz[0] = hz[rmi] // minus phase
				for rpi := 0; rpi < nhz; rpi++ {
					rphz[1] = hz[rpi] // plus phase

					ss.ResetTimePlots()
					ss.Init()
					for ri := 0; ri < ss.NReps; ri++ {
						Sp := float32(1)
						Rp := float32(1)
						for pi := 0; pi < nphs; pi++ {
							dur := phsdur[pi]
							shz := sphz[pi]
							rhz := rphz[pi]
							Sint := mat32.Exp(-1000.0 / float32(shz))
							Rint := mat32.Exp(-1000.0 / float32(rhz))
							for msec := 0; msec < dur; msec++ {
								Sp *= rand.Float32()
								prespike := false
								if Sp <= Sint {
									prespike = true
									Sp = 1
								}

								ge := float32(0.0)
								Rp *= rand.Float32()
								if Rp <= Rint {
									ge = ss.GeStim
									Rp = 1
								}
								if ss.RGClamp {
									ge = RGeStimForHz(float32(rhz))
								}

								ss.NeuronUpdt(msec, ge, 0, prespike)
								ss.LogDefault(0)
								if ss.StopNow {
									ss.Stopped()
									return
								}
							}
						}
						ss.LearnNow()
						ss.GraphRun(ss.ISISec, 0)
					}
					ss.GraphRun(ss.FinalSecs, 0)
					ss.LogPhaseDWt(ss.Log("PhaseDWtLog"), sphz, rphz)
					ss.Plot("PhaseDWtPlot").GoUpdate()
				}
			}
		}
	}
	ss.Stopped()
}
