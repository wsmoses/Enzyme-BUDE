SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
# .ONESHELL: # This doesn't work on all systems. Disabling for now.
.DELETE_ON_ERROR:

MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

ENZYME_PATH ?= /home/jhuckelheim/enzyme-sc22/Enzyme/enzyme/build/Enzyme/ClangEnzyme-15.so
CLANG_PATH ?= /home/jhuckelheim/enzyme-sc22/llvm-project/build/bin

# ENZYME_PATH ?= /home/wmoses/git/Enzyme/enzyme/buildomp/Enzyme/ClangEnzyme-14.so
# CLANG_PATH ?= /mnt/sabrent/wmoses/omp-llvm/build/bin

# ENZYME_PATH ?= /home/ubuntu/Enzyme/enzyme/buildD/Enzyme/ClangEnzyme-14.so
# CLANG_PATH ?= /home/ubuntu/llvm-project/buildD/bin


#ENZYME_PATH ?= /home/wmoses/git/Enzyme/enzyme/build13D/Enzyme/ClangEnzyme-13.so
#CLANG_PATH ?= /mnt/sabrent/wmoses/llvm13/buildall/bin

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

CFLAGS   = -std=c99 -Wall -mllvm -openmp-opt-max-iterations=1024 -mllvm -attributor-max-iterations=1024 -Ofast -ffast-math -$(ARCHFLAG)=$(ARCH) $(if $(AVX512),-mprefer-vector-width=512) -DWGSIZE=$(WGSIZE) -I $(OPENMP_PATH) -fno-exceptions -flegacy-pass-manager -Xclang -load -Xclang $(ENZYME_PATH) -Rpass=enzyme -mllvm -enzyme-print -mllvm -enzyme-max-int-offset=1024 -mllvm -enzyme-max-type-offset=1024  #-mllvm -enzyme-inline=1

LDFLAGS = -lm


# -------

.PHONY: all $(EXE) clean

all: omp-single-forward.exe ser-single-forward.exe omp-single-gradient.exe ser-single-gradient.exe ompOpt-single-forward.exe ompOpt-single-gradient.exe

verify: omp-single-forward-verify.exe omp-single-gradient-verify.exe ompOpt-single-gradient-verify.exe
	./ompOpt-single-gradient-verify.exe | grep Checksum
	./omp-single-gradient-verify.exe | grep Checksum
	./omp-single-forward-verify.exe | grep Checksum

# all: omp-single-forward.exe omp-single-gradient.exe

ser-single-forward.exe: $(SOURCES2.0)
	time $(CC) bude.c $(CFLAGS) -lm -o $@ -DFORWARD=1

omp-single-forward.exe: $(SOURCES2.0)
	time $(CC) -fopenmp bude.c $(CFLAGS) -lm -o $@ $(OPENMP_LIB) -DFORWARD=1 -mllvm -enzyme-omp-opt=0 

omp-single-forward-verify.exe: $(SOURCES2.0)
	time $(CC) -fopenmp bude.c $(CFLAGS) -lm -o $@ $(OPENMP_LIB) -DFORWARD=1 -mllvm -enzyme-omp-opt=0 -DVERIFY=1

ser-single-gradient.exe: $(SOURCES2.0)
	time $(CC) bude.c $(CFLAGS) -lm -o $@

omp-single-gradient-verify.exe: $(SOURCES2.0)
	$(CC) -fopenmp bude.c $(CFLAGS) -lm -o $@ $(OPENMP_LIB) -mllvm -enzyme-omp-opt=0 -g -DVERIFY=1

omp-single-gradient.exe: $(SOURCES2.0)
	time $(CC) -fopenmp bude.c $(CFLAGS) -lm -o $@ $(OPENMP_LIB) -mllvm -enzyme-omp-opt=0 

ompOpt-single-forward.exe: $(SOURCES2.0)
	time $(CC) -fopenmp bude.c $(CFLAGS) -lm -o $@ $(OPENMP_LIB) -DFORWARD=1 -mllvm -enzyme-omp-opt=1

ompOpt-single-gradient.exe: $(SOURCES2.0)
	time $(CC) -fopenmp bude.c $(CFLAGS) -lm -o $@ $(OPENMP_LIB) -mllvm -enzyme-omp-opt=1

ompOpt-single-gradient-verify.exe: $(SOURCES2.0)
	time $(CC) -fopenmp bude.c $(CFLAGS) -lm -o $@ $(OPENMP_LIB) -mllvm -enzyme-omp-opt=1 -DVERIFY=1

clean:
	rm -f *.exe

