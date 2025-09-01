JL = julia --project

default: init update

init:
	$(JL) -e 'using Pkg; Pkg.instantiate()'

update:
	$(JL) -e 'using Pkg; Pkg.update(); Pkg.precompile()'