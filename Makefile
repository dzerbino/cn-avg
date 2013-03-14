PREFIX = /inside/home/dzerbino/.local
export BIN_DIR = ${PWD}/bin
export LIB_DIR = ${PWD}/lib
export SONLIB_DIR = /inside/home/dzerbino/utils/sonLib
export CACTUS_DIR = /inside/home/dzerbino/utils/cactus

default:
	cd src && make
	cd bin && chmod 755 *

install:
	cp ${BIN_DIR}/* ${PREFIX}/bin
	cp ${LIB_DIR}/* ${PREFIX}/lib

clean:
	rm ${BIN_DIR}/*
	rm ${LIB_DIR}/*

doc: cnavg
	epydoc --graph all --no-private --html --output doc/html cnavg
