SHELL = bash
COMP = g++
CCOMP = gcc
PREFIX ?= /usr/local
CXXFLAGS = -std=c++11 -fPIC -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2)
CXXIFLAGS = -I$(PREFIX)/include -Iinclude
CIFLAGS = -I$(PREFIX)/include -Iinclude
LFLAGS = -L$(PREFIX)/lib -Llib

ifeq ($(findstring align_pipelines, ${CONDA_PREFIX}), align_pipelines)
    CXXIFLAGS += -I${CONDA_PREFIX}/include
    CIFLAGS += -I${CONDA_PREFIX}/include
    LFLAGS += -L${CONDA_PREFIX}/lib
endif

BC_LENX2 = 32
KX2 = 16
DEPS = lib/libhtswrapper.a
DEPS2 = -lz -lhts -lpthread

atac_fq_preprocess: src/atac_fq_preprocess.cpp $(DEPS)
	$(COMP) $(CXXFLAGS) $(CXXIFLAGS) src/atac_fq_preprocess.cpp $(LFLAGS) $(DEPS) -o atac_fq_preprocess $(DEPS2)

lib/libhtswrapper.a:
	cd dependencies/htswrapper && $(MAKE) PREFIX=../..
	cd dependencies/htswrapper && $(MAKE) install PREFIX=../..

clean: clean_deps
	rm -f atac_fq_preprocess
	rm lib/libhtswrapper.a

clean_deps:
	cd dependencies/htswrapper && $(MAKE) clean || true

