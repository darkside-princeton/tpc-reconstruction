TARGET = waveTPC
PROGNAME = waveTPC

IDIR = inc
SDIR = src
VPATH = src

DEPS = $(wildcard $(IDIR)/*.h)
DEPS += $(wildcard $(SDIR)/*.cc)

_SRCS = $(wildcard $(SDIR)/*.cc)
_OBJS = $(_SRCS:.cc=.o)
OBJS = $(patsubst $(IDIR)/%,$(SDIR)/%,$(_OBJS))
TOBJ = $(TARGET).o

LDFLAGS = -g
LIBS = $(shell root-config --glibs)

CC = g++
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
CFLAGS = -O3 -c -g -Wall -Wno-write-strings -I$(IDIR)
COPTS = -fPIC -DLINUX -Wall \
        $(shell root-config --cflags)

$(TARGET): $(OBJS) $(TOBJ) $(DEPS)
	@echo [g++] Linking...
	@$(CC) $(TOBJ) $(OBJS) -o $(PROGNAME) $(LDFLAGS) $(LIBS) -L$(IFDHC_LIB) -lifdh -lSpectrum
	@echo [DONE]

$(SDIR)/%.o: src/%.cc $(DEPS)
	@echo [g++] Compiling class $<
	@$(CC) $(CFLAGS) $(COPTS) $< -o $@ $(LDFLAGS) $(LIBS)

$(TOBJ) : $(TARGET).cc $(DEPS)
	@echo [g++] Compiling $(TARGET).cc
	@$(CC) $(CFLAGS) $(COPTS) $< -o $@ $(LDFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -rf $(SDIR)/*.o *.o $(TOBJ)

clobber:
	rm -rf $(SDIR)/*.*~ $(IDIR)/*.*~ *.*~  
