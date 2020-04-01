package main

import (
	"os"
	"path/filepath"

	clicmdflags "github.com/codemodify/systemkit-clicmdflags"
)

// AppRootCmdFlags -
type AppRootCmdFlags struct {
	Method            string `flagName:"method"            flagRequired:"true" flagDescription:"Estimation approach, \"iid\" or \"non-iid\""`
	File              string `flagName:"file"              flagRequired:"true" flagDescription:"Binary file with at least 1 million entries (samples)"`
	BitsPerSample     int    `flagName:"bitsPerSample"     flagRequired:"true" flagDescription:"Must be between 1-8, inclusive"`
	InitialEE         bool   `flagName:"initialEE"         flagDefault:"true"  flagDescription:"For initial entropy estimate (Section 3.1.3)"`
	ConditionedEE     bool   `flagName:"conditionedEE"     flagDefault:"false" flagDescription:"For conditioned sequential dataset entropy estimate (Section 3.1.5.2)"`
	TestAllBits       bool   `flagName:"testAllBits"       flagDefault:"true"  flagDescription:"Tests all bits in bitstring"`
	TruncateBitstring bool   `flagName:"truncateBitstring" flagDefault:"false" flagDescription:"Truncates bitstring to 1000000 bits"`
	ReadSubstring     string `flagName:"readSubstring"     flagDefault:"0-10"  flagDescription:"Read the <index> substring of length <samples>"`
	JSON              bool   `flagName:"json"              flagDefault:"false" flagDescription:"Enables JSON output"`
	Verbose           bool   `flagName:"verbose"           flagDefault:"false" flagDescription:"Enables verbose output"`
}

var appRootCmd = &clicmdflags.Command{
	Name:        filepath.Base(os.Args[0]),
	Description: "Entropy Assessment",
	Examples: []string{
		filepath.Base(os.Args[0]) + " -method iid     -file FILE_TO_TEST",
		filepath.Base(os.Args[0]) + " -method non-iid -file FILE_TO_TEST",
	},
	Flags: AppRootCmdFlags{},
}

// Execute - this is a convenience call
func Execute() error {
	return appRootCmd.Execute()
}
