#Compiler
CXX 	 := g++

#Executable program
TARGET 	 :=	DarkAxion

#Source, build and executable directories
SRCDIR    := SourceFiles
BUILDDIR  := BuildFiles
TARGETDIR := BinFiles
INPUTDIR  := InputFiles
OUTPUTDIR := OutputFiles
OUTPUTDIR += OutputFiles/PNGPlots
OUTPUTDIR += OutputFiles/PNGPlots/Variables
OUTPUTDIR += OutputFiles/PNGPlots/VariablesComparisonForDifferentProductions
OUTPUTDIR += OutputFiles/PNGPlots/VariablesComparisonForDifferentProductionsHW
OUTPUTDIR += OutputFiles/PNGPlots/VariablesComparisonForDifferentProductionsHW/OneStripping
OUTPUTDIR += OutputFiles/PNGPlots/VariablesComparisonForDifferentProductionsHW/TwoStripping
OUTPUTDIR += OutputFiles/PNGPlots/Angles
OUTPUTDIR += OutputFiles/PNGPlots/Angles/Quadrati
OUTPUTDIR += OutputFiles/PNGPlots/Discriminants
OUTPUTDIR += OutputFiles/PNGPlots/Efficiencies
OUTPUTDIR += OutputFiles/PNGPlots/Efficiencies/Q_squared
OUTPUTDIR += OutputFiles/PNGPlots/Efficiencies/Missing_mass_squared
OUTPUTDIR += OutputFiles/PNGPlots/Efficiencies/Missing_pt
OUTPUTDIR += OutputFiles/PNGPlots/Efficiencies/ctau
OUTPUTDIR += OutputFiles/PNGPlots/Efficiencies/Q_squared_quadrati
OUTPUTDIR += OutputFiles/PNGPlots/Efficiencies/Missing_mass_squared_quadrati
OUTPUTDIR += OutputFiles/PNGPlots/Efficiencies/Missing_pt_quadrati
OUTPUTDIR += OutputFiles/PNGPlots/Efficiencies/ctau_quadrati
OUTPUTDIR += OutputFiles/PNGPlots/FitResults
OUTPUTDIR += OutputFiles/TMVAResults/ClosureTest
OUTPUTDIR += OutputFiles/TMVAResults/ClosureTest/Normalized
OUTPUTDIR += OutputFiles/TMVAResults/ClosureTest/Normalized/Pass_NNFail_Ratio

#Source and build file extension
SRCEXT := cpp
OBJEXT := o
DEPEXT := d

#Flags
CXXFLAGS 	 := -O3 -g -Wall $(shell root-config --cflags)

#Includes
INCROOT 	 := -I$(shell root-config --incdir)
INCROOT 	 += -I/home/pasquale/Blue-2.4.0
INCROOT 	 += -I/home/root/include/Math
INCROOT 	 += -I/home/root/include/RooStats
INCLUDES     := $(INCROOT)

#Libraries that has to be linked
LIBSROOT     := -L/home/pasqualeandreola/Blue-2.4.0
LIBSROOT 	 += $(shell root-config --libs) -lRooFitCore -lRooFit -lTMVA -lTMVAGui -lXMLIO -lMathMore -lSPlot -lRooStats 
LIBS 	 	 := $(LIBSROOT)

#Sources and objects research
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))	

#Default make execution
all: directories $(TARGET)
	 @./BinFiles/DarkAxion

#Remake
remake: cleaner all

#Making of directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)
	@mkdir -p $(INPUTDIR)
	@mkdir -p $(OUTPUTDIR)

#Clean Build
clean:
	@$(RM) -rf $(BUILDDIR)

#Clean everything
cleaner: clean
	@$(RM) -rf $(TARGETDIR)

#Documenting
doxy:
	doxygen Documentation/doxy_config

#Testing
test:
	@echo 
	@echo To make tests, open the TestFiles directory and run the command "make"\n
	@echo
	
#Getting the dipendencies
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Setting up the linker
$(TARGET): $(OBJECTS)
	$(CXX) -o $(TARGETDIR)/$(TARGET) $^ $(LIBS)

#Compiling
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp
	@rm -rf $(BUILDDIR)/$*.$(DEPEXT)

.PHONY: all remake clean cleaner doxy
