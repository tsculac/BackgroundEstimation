#makefile

CXXFLAGS = -g -I. -m64 $(shell root-config --cflags) -I include
LDFLAGS = $(shell root-config --libs) -lm -lGenVector
CXX = g++

#LIBS = /afs/cern.ch/work/m/mkovac/CMS/RUN2/HZZ4l_Plotter/CMSSW_8_0_21/src/ext/
#EXTLIBS = -L$(LIBS) cConstants_cc.so FinalStates_cc.so bitops_cc.so

EXTLIBS = ./ext/cConstants_cc.so ./ext/FinalStates_cc.so ./ext/bitops_cc.so

VPATH = ./src/

SRCPP = run.cpp\
        OSmethod.cpp\
		  Tree.cpp\
		  Settings.cpp\
		  Category.cpp
         
OBJCPP = $(patsubst %.cpp,obj/%.o,$(SRCPP))

all : run

obj/%.o : %.cpp
	@echo ">> compiling $*"
	@mkdir -p obj/
	@$(CXX) -c $< ${CXXFLAGS} -o $@

run : $(OBJCPP)
	@echo ">> linking..."
	@$(CXX) $^ $(EXTLIBS) ${LDFLAGS} ${CXXFLAGS}  -o $@

clean:
	@echo ">> cleaning objects and executable"
	@rm  -f obj/*.o
	@rm -f run
	
uninstall:
	@echo ">> Uninstalling plotter"
	@rm  -f obj/*.o
	@rm  -f ext/*.so ext/*.d ext/*.pcm
	@rm -f run
