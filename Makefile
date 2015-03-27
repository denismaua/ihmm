MAKE = make
SRCDIR = src
EXEDIR = bin
OBJDIR = src

.PHONY: build clean dist-clean

all: build

build:
	$(MAKE) -C $(SRCDIR)

filter:
	$(MAKE) -C $(SRCDIR) filter

likelihood:
	$(MAKE) -C $(SRCDIR) likelihood

viterbi:
	$(MAKE) -C $(SRCDIR) viterbi

filter2:
	$(MAKE) -C $(SRCDIR) filter2


clean:
	$(MAKE) -C $(SRCDIR) clean

dist-clean:
	$(MAKE) -C $(SRCDIR) dist-clean
