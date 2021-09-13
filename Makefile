CODEINFO = -DCODEFLAG="\"official version\"" \
	   -DCODEMAJORVER="\"1\"" \
	   -DCODEMINORVER="\"0\"" \
	   -DCODEUSER="\"`whoami | sed 's/\\\\/\\\\\\\\/g'`\"" \
	   -DCODEMACHINE="\"`hostname`\"" \
	   -DCODECOMPILER="\"`$(CXX) --version | head -n 1`"\" \
           -DCODECOMPFLAG="\"$(CXXFLAG)"\" \
	   -DCODELINUX \
       -DEIGEN_DONT_VECTORIZE

CXX     = g++
CXXFLAG =-O3 -Os -fno-exceptions --std=c++11 -fopenmp #  -static-libstdc++ -static-libgcc -static -mavx -msse4 #-lboost_system -Wall -Wno-unknown-pragmas -Wno-sign-compare -Wno-unused -g
CC      = gcc
CCFLAG  =-O3 -fopenmp -Wall # -g

INCDIR  =./include
LIBDIR  =./lib
SRCS    =./src
BINS    =./bin
INCS    =-I$(INCDIR) -I$(SRCS) -I/home/zhan055/lib 
LIBS    =-L$(LIBDIR) #-lpthread -llapack -lblas
OBJS    =./obj
DEPS    =./depend

TARGET   = $(BINS)/nwpesse
SOURCES := $(wildcard $(SRCS)/*.cpp)
OBJECTS := $(patsubst $(SRCS)%.cpp,$(OBJS)%.o,$(SOURCES))    

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAG) $(INCS) $(LIBS) -o $@ $^


sinclude $(patsubst $(SRCS)%.cpp,$(DEPS)%.d,$(SOURCES))

$(OBJS)/%.o: $(SRCS)/%.cpp
	$(CXX) $(CXXFLAG) $(CODEINFO) $(INCS) $(LIBS) -c -o $@ $< 


$(DEPS)/%.d: $(SRCS)/%.cpp
	@set -e; \
	rm -f $@; \
	$(CXX) -MM $(CXXFLAG) $(CODEINFO) $(INCS) $< > $@.$$$$;\
	sed 's/\([a-zA-Z0-9]*\)\.o[ :]/$(patsubst ./%,%,$(OBJS)\/$(shell v1=$@; v2=$${v1%/*}; v3=$${v2#*/}; echo $${v3} | sed 's/\//\\\//g'))\/\1.o $(patsubst ./%,%,$(DEPS)\/$(shell v1=$@; v2=$${v1%/*}; v3=$${v2#*/}; echo $${v3} | sed 's/\//\\\//g'))\/\1.d:/g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

.PHONY: clean
clean:
	rm -rf $(OBJECTS) $(TARGET) $(DEPS)/*.d*

