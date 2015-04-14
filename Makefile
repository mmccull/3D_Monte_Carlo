

CXX = g++
LAPACK = /Users/mmccull/lib
OPTS = -O3 -ftree-vectorize

mc: liquid_argon_mc.c stringlib.c
	$(CC) -c  liquid_argon_mc.c stringlib.c $(OPTS) 
	$(CC) liquid_argon_mc.o stringlib.o $(OPTS) -o liquid_argon_mc.x
#	cp sasa_force.x ../3app

