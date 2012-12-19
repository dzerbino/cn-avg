PREFIX = /inside/home/graim/.local

default:
	cd src && make
	cd bin && chmod 755 *

install:
	cp bin/* ${PREFIX}/bin

clean:
	rm bin/*

doc: cnavg
	epydoc --graph all --no-private --html --output doc/html cnavg
