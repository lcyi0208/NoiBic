PROGS = noibic
SRCS  = main.cpp struct.cpp get_options.cpp data_processing.cpp cluster.cpp cluster_expand.cpp LCS.cpp result_output.cpp seed_cal.cpp seed_generation.cpp
OBJS  = $(SRCS:.cpp=.o)
CC    = g++
BOOST_INCLUDE ?= /usr/include

LDFLAGS  = -lm -pthread
CPPFLAGS = -std=c++0x -O3 -g -I $(BOOST_INCLUDE) -pthread


VERSION ?= 1.0.1

GIT_COMMIT := $(shell git rev-parse --short HEAD 2>/dev/null)
GIT_TAG    := $(shell git describe --tags --abbrev=0 2>/dev/null)
GIT_DIRTY  := $(shell test -n "$$(git status --porcelain 2>/dev/null)" && echo "-dirty")
BUILD_DATE := $(shell date -u '+%Y-%m-%dT%H:%M:%SZ')
HOST       := $(shell uname -s)-$(shell uname -m)
VCS_VER    := $(if $(GIT_TAG),$(GIT_TAG),$(GIT_COMMIT))
FULL_VERSION := $(VERSION)$(if $(VCS_VER),+$(VCS_VER)$(GIT_DIRTY),)


CPPFLAGS += \
  -DAPP_NAME=\"$(PROGS)\" \
  -DAPP_VERSION=\"$(FULL_VERSION)\" \
  -DAPP_BUILD_DATE=\"$(BUILD_DATE)\" \
  -DAPP_HOST=\"$(HOST)\"


.PHONY: all clean version
all: $(PROGS)

$(PROGS): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $< -o $@

version:
	@echo "$(PROGS) version: $(FULL_VERSION)"
	@echo "built at: $(BUILD_DATE) on $(HOST)"

clean:
	rm -f $(PROGS) *.o
// Initial commit (no-op)
