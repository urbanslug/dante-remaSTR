SHELL := bash
.ONESHELL:
.SHELLFLAGS := -euo pipefail -c
.DELETE_ON_ERROR:

.PHONY: linux_build
linux_build:
	mkdir -p build
	cd build
	pyinstaller -F ../dante_remastr_standalone.py
	cd dist
	cp -r ../../dante_remastr_standalone_templates .
	cp -r ../../includes .
