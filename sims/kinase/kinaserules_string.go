// Code generated by "stringer -type=KinaseRules"; DO NOT EDIT.

package main

import (
	"errors"
	"strconv"
)

var _ = errors.New("dummy error")

func _() {
	// An "invalid array index" compiler error signifies that the constant values have changed.
	// Re-run the stringer command to generate them again.
	var x [1]struct{}
	_ = x[NeurSpkCa-0]
	_ = x[SynSpkCa-1]
	_ = x[SynSpkNMDA-2]
	_ = x[SynNMDACa-3]
	_ = x[KinaseRulesN-4]
}

const _KinaseRules_name = "NeurSpkCaSynSpkCaSynSpkNMDASynNMDACaKinaseRulesN"

var _KinaseRules_index = [...]uint8{0, 9, 17, 27, 36, 48}

func (i KinaseRules) String() string {
	if i < 0 || i >= KinaseRules(len(_KinaseRules_index)-1) {
		return "KinaseRules(" + strconv.FormatInt(int64(i), 10) + ")"
	}
	return _KinaseRules_name[_KinaseRules_index[i]:_KinaseRules_index[i+1]]
}

func (i *KinaseRules) FromString(s string) error {
	for j := 0; j < len(_KinaseRules_index)-1; j++ {
		if s == _KinaseRules_name[_KinaseRules_index[j]:_KinaseRules_index[j+1]] {
			*i = KinaseRules(j)
			return nil
		}
	}
	return errors.New("String: " + s + " is not a valid option for type: KinaseRules")
}
