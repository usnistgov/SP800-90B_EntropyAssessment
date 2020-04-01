package main

import (
	"fmt"
	"os"
	"path/filepath"

	clicmdflags "github.com/codemodify/systemkit-clicmdflags"
)

func init() {
	appRootCmd.AddCommand(&clicmdflags.Command{
		Name:        "version",
		Description: "Displays version",
		Examples: []string{
			filepath.Base(os.Args[0]) + " version",
		},
		Handler: func(command *clicmdflags.Command) {
			fmt.Println("v1.0")
		},
	})
}
