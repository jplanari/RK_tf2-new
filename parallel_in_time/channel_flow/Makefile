# ------------------------------------------------------------------------------
# Tools
# ------------------------------------------------------------------------------
MAKEMOD=tf-makemod

# ------------------------------------------------------------------------------
# Flags
# ------------------------------------------------------------------------------
MAKEFLAGS += --no-builtin-rules

# ------------------------------------------------------------------------------
# Auto computed flags
# ------------------------------------------------------------------------------
CXXFLAGS=-isystem /usr/lib/x86_64-linux-gnu/openmpi/include
CXXFLAGS+=-Wall -Wextra -Wconversion -Wsign-conversion -Wfloat-equal
export CXXFLAGS

################################################################################
# Automatic part
################################################################################

MODS_DIR=./specific_mods
MODS=$(shell find ./$(MODS_DIR) -maxdepth 1 -name "*.cpp")
SOBJS=$(patsubst %.cpp,%.so,$(MODS))

# ******************************************************************************
# Targets
# ******************************************************************************

.DEFAULT_GOAL := all

.PHONY: all debug
all: $(SOBJS)

# This enables 'make debug' to be equal to 'make' but with debugging info on.
debug: DEBUG=-d
debug: $(SOBJS)

%.so: %.cpp $(TFLIB_DEP)
	$(MAKEMOD) $< $(DEBUG)

# ******************************************************************************
# Unconditional targets
# ******************************************************************************
.PHONY: clean install uninstall
clean:
	@rm -f $(SOBJS)
