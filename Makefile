# ===== Makefile for call_features =====
# Usage examples:
#   make               # release build, optimised
#   make DEBUG=1       # debug build
#   make OPENMP=1      # enable OpenMP (add DEBUG=1 if needed)
#   make OPENMP=1 THREADS=16   # build with OpenMP, will run with 16 threads by default
#   OMP_NUM_THREADS still overrides at run-time
#   make clean         # remove binary

# ---------------------------------------------------------------------------
# Toolchain
CC      ?= gcc

# ---------------------------------------------------------------------------
# Project files / output
SRCS    := src/flex_assign.c
BINARY  := flex_demux_mtx

# ---------------------------------------------------------------------------
# Common flags
WARN    := -Wall -Wextra -pedantic
STD     := -std=c11
OPT     := -O3 -march=native
DBG     := -g -O0 -DDEBUG     # -DDEBUG lets you gate extra debug code if desired

# ---------------------------------------------------------------------------
# Build-type selection
ifeq ($(DEBUG),1)
    CFLAGS := $(WARN) $(STD) $(DBG)
else
    CFLAGS := $(WARN) $(STD) $(OPT)
endif

# ---------------------------------------------------------------------------
# Optional OpenMP support
ifeq ($(OPENMP),1)
    CFLAGS  += -fopenmp -DUSE_OPENMP
    LDLIBS  += -fopenmp
    # allow a default thread count at link-time (can be overridden later)
    ifneq ($(THREADS),)
        CFLAGS += -DDEFAULT_OMP_THREADS=$(THREADS)
    endif
endif

# ---------------------------------------------------------------------------
# Libraries
LDLIBS  += -lm

# ---------------------------------------------------------------------------
# Executable output directory
BINDIR := bin

# Executable 1 – Flex Demux (fast, single-file)
BIN1   := flex_demux_mtx
SRC1   := src/flex_assign.c

# Executable 2 – Call Features (streaming + optional EM + OpenMP)
BIN2   := call_features
SRC2   := src/call_features.c src/per_guide_em.c src/compute_Mmin_from_cells.c
INCS   := -Iincludes          # header path for per_guide_em.h and compute_Mmin_from_cells.h

# Paths of final binaries
OUT1 := $(BINDIR)/$(BIN1)
OUT2 := $(BINDIR)/$(BIN2)

# ---------------------------------------------------------------------------
# Common build rule macro
define build_template
$(1): $(2)
	@mkdir -p $(BINDIR)
	$$(CC) $$(CFLAGS) $$(INCS) $$^ -o $$@ $$(LDLIBS)
endef

# Instantiate the two rules
$(eval $(call build_template,$(OUT1),$(SRC1)))
$(eval $(call build_template,$(OUT2),$(SRC2)))

# ---------------------------------------------------------------------------
.PHONY: all clean
all: $(OUT1) $(OUT2)

clean:
	rm -f $(OUT1) $(OUT2)
# ---------------------------------------------------------------------------
