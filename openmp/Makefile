SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
# .ONESHELL: # This doesn't work on all systems. Disabling for now.
.DELETE_ON_ERROR:

MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

ENZYME_PATH ?= /home/wmoses/git/Enzyme/enzyme/build13D/Enzyme/ClangEnzyme-13.so
CLANG_PATH ?= /mnt/sabrent/wmoses/llvm13/buildall/bin

OPENMP_PATH ?= $(CLANG_PATH)/../projects/openmp/runtime/src
MPI_PATH ?= /usr/lib/x86_64-linux-gnu/openmpi/include
OPENMP_LIB ?= $(CLANG_PATH)/../lib/libomp.so

# -------

WGSIZE = 64
ARCH   = native
AVX512 =
EXTRA_FLAGS =

# -------

MACHINE = $(shell uname -m)

ifeq ($(MACHINE), x86_64)
ARCHFLAG = march

# On Skylake and Cascade Lake, 256-bit vectors are used by default, but 512-bit is more beneficial
AVX512   =  $(if $(filter skylake% cascadelake% %avx512, $(ARCH)),yes,)

# Intel platforms benefit more from unrolling, specially Skylake and later
WGSIZE = 256

else
# The target CPU is specificed differently on x86 and on aarch64
# https://community.arm.com/developer/tools-software/tools/b/tools-software-ides-blog/posts/compiler-flags-across-architectures-march-mtune-and-mcpu
ARCHFLAG = mcpu
endif

# -------

CC = $(CLANG_PATH)/clang

CFLAGS   = -std=c99 -Wall -Ofast -ffast-math -$(ARCHFLAG)=$(ARCH) $(if $(AVX512),-mprefer-vector-width=512) -DWGSIZE=$(WGSIZE) -I $(OPENMP_PATH) -fno-exceptions -fno-experimental-new-pass-manager -Xclang -load -Xclang $(ENZYME_PATH) -Rpass=enzyme -mllvm -enzyme-print -mllvm -enzyme-max-int-offset=1024 -mllvm -enzyme-max-type-offset=1024  #-mllvm -enzyme-inline=1

LDFLAGS = -lm


# -------

.PHONY: all $(EXE) clean

all: omp-single-forward.exe ser-single-forward.exe omp-single-gradient.exe ser-single-gradient.exe

# all: omp-single-forward.exe omp-single-gradient.exe

ser-single-forward.exe: $(SOURCES2.0)
	time $(CC) bude.c $(CFLAGS) -lm -o $@ -DFORWARD=1

omp-single-forward.exe: $(SOURCES2.0)
	time $(CC) -fopenmp bude.c $(CFLAGS) -lm -o $@ $(OPENMP_LIB) -DFORWARD=1

ser-single-gradient.exe: $(SOURCES2.0)
	time $(CC) bude.c $(CFLAGS) -lm -o $@

omp-single-gradient.exe: $(SOURCES2.0)
	time $(CC) -fopenmp bude.c $(CFLAGS) -lm -o $@ $(OPENMP_LIB)

clean:
	rm -f *.exe

