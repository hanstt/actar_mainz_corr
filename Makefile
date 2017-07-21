all: analyze corr recon

analyze: analyze.C Makefile
	eg++ -std=c++11 -Wall `root-config --cflags` -o $@ $< `root-config --libs` `pkg-config freetype2 --libs`

corr: corr.C Makefile
	eg++ -std=c++11 -Wall `root-config --cflags` -o $@ $< `root-config --libs` `pkg-config freetype2 --libs`

recon: recon.C Makefile
	eg++ -std=c++11 -Wall `root-config --cflags` -o $@ $< `root-config --libs` `pkg-config freetype2 --libs`
