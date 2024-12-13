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
	cp -r ../../templates .
	cp -r ../../includes .

# https://stackoverflow.com/questions/2950971/packaging-a-python-script-on-linux-into-a-windows-executable
