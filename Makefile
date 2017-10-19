#makefile

CXXFLAGS = -g -I. -m64 $(shell root-config --cflags) -I include
LDFLAGS = $(shell root-config --libs) -lm -lGenVector
CXX = g++

EXTLIBS = ./ext/cConstants_cc.so ./ext/FinalStates_cc.so ./ext/bitops_cc.so ./ext/Discriminants_cc.so

VPATH = ./src/ ./include/

SRCPP_OS = run_OS.cpp\
        OSmethod.cpp\
		  Tree.cpp\
		  Settings.cpp\
		  Category.cpp\
		  FakeRates.cpp\
		  Plots.cpp\
		  CMS_lumi.cpp

SRCPP_SS = run_SS.cpp\
        SSmethod.cpp\
		  Tree.cpp\
		  Settings.cpp\
		  Category.cpp\
		  FakeRates.cpp\
		  Plots.cpp\
		  CMS_lumi.cpp

INCLUDES = OSmethod.h\
	   SSmethod.h\
	   Tree.h\
		  Settings.h\
		  Category.h\
		  FakeRates.h\
		  Plots.h\
		  CMS_lumi.h\

         
OBJCPP_OS = $(patsubst %.cpp,obj/%.o,$(SRCPP_OS))
OBJCPP_SS = $(patsubst %.cpp,obj/%.o,$(SRCPP_SS))

all : run_OS run_SS

obj/%.o: %.cpp $(INCLUDES)
	@echo ">> compiling $*"
	@mkdir -p obj/
	@$(CXX) -c $< ${CXXFLAGS} -o $@

run_OS : $(OBJCPP_OS)
	@echo ">> linking..."
	@$(CXX) $^ $(EXTLIBS) ${LDFLAGS} ${CXXFLAGS}  -o $@

run_SS : $(OBJCPP_SS)
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
