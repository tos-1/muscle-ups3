all: build

#SRC := src/halomodel.c

.PHONY: build clean

build:
	python setup.py build_ext --inplace

clean:
	- rm -f ext/*
