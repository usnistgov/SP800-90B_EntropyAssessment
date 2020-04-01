package main

import (
	"log"
)

func main() {
	if err := Execute(); err != nil {
		log.Fatal(err)
	}
}
