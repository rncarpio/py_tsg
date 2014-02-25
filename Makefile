
DEBUG = 0
VPATH = TasmanianSparseGrids

ifeq ($(OS),Windows_NT)
	COMPILER = msvc
else
	COMPILER = gcc
endif

#############################
# include directories
#############################

ifeq ($(COMPILER), gcc)		
	BOOST_DIR = $(HOME)/boost_1_54_0
	BOOST_INC_DIR = $(BOOST_DIR)
	PYTHON_DIR = $(HOME)/local
	PYTHON_INC_DIRS = $(HOME)/local/include/python2.7 $(HOME)/local/lib/python2.7/site-packages/numpy/core/include
else	
	BOOST_DIR = c:/boost/boost_1_55_0
	BOOST_INC_DIR = $(BOOST_DIR)
	PYTHON_DIR = c:/Python27
	PYTHON_INC_DIRS = $(PYTHON_DIR)/include $(PYTHON_DIR)/lib/site-packages/numpy/core/include
	PYTHON_LIB_DIRS = $(PYTHON_DIR)/libs
endif
PYUBLAS_INC_DIR = .
INCLUDES = -Isrc -I$(BOOST_INC_DIR) -I$(PYUBLAS_INC_DIR) \
	$(foreach dir, $(PYTHON_INC_DIRS), -I$(dir))

#############################
# library directories
#############################
	
BOOST_LIB_DIR = $(BOOST_DIR)/lib
ifeq ($(COMPILER), gcc)	
	PYTHON_LIB = python2.7
	BOOST_PYTHON_LIB = boost_python
	LIB_PATHS = $(HOME)/local/lib $(HOME)/local/src/boost_1_51_0/stage/lib
	LIBS = -l$(BOOST_PYTHON_LIB) -l$(PYTHON_LIB)
	LIB_DIRS = $(foreach dir, $(LIB_PATHS), -L$(dir))
else	
	PYTHON_LIB = python27
	BOOST_PYTHON_LIB = boost_python-vc100-mt-1_55
	LIB_PATHS = $(BOOST_LIB_DIR) $(PYTHON_LIB_DIRS)
	LIB_NAMES = $(BOOST_PYTHON_LIB) $(PYTHON_LIB)
	LIB_DIRS = $(foreach dir, $(LIB_PATHS), /LIBPATH:$(dir))
	LIBS = $(BOOST_LIB_DIR)/$(BOOST_PYTHON_LIB).lib $(PYTHON_LIB).lib
endif

BUILD_DIR = build
ifeq ($(COMPILER), msvc)	
	EXE_PREFIX = .exe
	DL_SUFFIX = .dll
else
	EXE_PREFIX = 
	DL_SUFFIX = .so
endif

#############################
# compile/link flags
#############################

ifeq ($(COMPILER), gcc)			
	CC = gcc
	CXX = g++ 
	LINK = g++	
	LINKFLAGS = -shared $(LINKFLAGS_EXE) -fopenmp
	ifeq ($(DEBUG), 0)
		CXXFLAGS = -O3 -ffast-math -mtune=native -fPIC -fopenmp
	else
		CXXFLAGS = -g -ffast-math -mtune=native -fPIC -fopenmp
	endif
else ifeq ($(COMPILER), msvc)
	CC = cl.exe
	CXX = cl.exe
#--compiler-options /MD
	LINK = link.exe
	ifeq ($(DEBUG), 0)
		CXXFLAGS = /nologo /O2 /MD /W3 /GS- /Zi /EHsc /Zm1000 /DNDEBUG /WL /fp:strict /DFP_STRICT
		LINKFLAGS_EXE = /nologo /INCREMENTAL:NO /DEBUG
	else
		CXXFLAGS = /nologo /Od /MD /W3 /GS- /Zi /EHsc /Zm1000 /WL /fp:strict /DFP_STRICT
		LINKFLAGS_EXE = /nologo /INCREMENTAL:NO /DEBUG
	endif
	LINKFLAGS = /DLL $(LINKFLAGS_EXE)
endif #msvc

#########################
## generic rules
#########################

ifeq ($(COMPILER), gcc)

$(BUILD_DIR)/%.obj: %.cpp %.hpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

$(BUILD_DIR)/%.obj: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

_%$(DL_SUFFIX):	$(BUILD_DIR)/%.o
	$(LINK) -o $@ $< $(LINKFLAGS) $(LIB_DIRS) $(LIBS)

$(BUILD_DIR)/%:	$(BUILD_DIR)/%.o
	$(LINK) $(LINKFLAGS_EXE) $(LIB_DIRS) $(LIBS) $< -o $@

else ifeq ($(COMPILER), msvc)

$(BUILD_DIR)/%.obj: %.cpp %.hpp
	$(CXX) /c $(CXXFLAGS) $(INCLUDES) /Tp$< /Fo$@
	
$(BUILD_DIR)/%.obj: %.cpp
	$(CXX) /c $(CXXFLAGS) $(INCLUDES) /Tp$< /Fo$@

_%$(DL_SUFFIX): $(BUILD_DIR)/%.obj 
	$(LINK) $(LINKFLAGS) $(LIB_DIRS) $(LIBS) $< /OUT:$@ \
	IMPLIB:$(BUILD_DIR)/_$*.lib MANIFESTFILE:$(BUILD_DIR)/_$*.$(DL_SUFFIX).manifest

$(BUILD_DIR)/%.exe: $(BUILD_DIR)/%.obj
	$(LINK) $(LINKFLAGS_EXE) $(LIB_DIRS) $(LIBS) $< /OUT:$@
	
endif

#########################
## targets
#########################

LHEADERS = tsgIndexSet.hpp tsgHelperFunctions.hpp tsgEnumerate.hpp \
		tsgBase1DRule.hpp tsgRuleClenshawCurtis.hpp tsgRuleChebyshev.hpp tsgRuleGaussLegendre.hpp tsgRulePieceWiseLocal.hpp \
		tsgRuleChebyshevNestedTwoPoint.hpp tsgRuleGaussChebyshev1.hpp tsgRuleGaussChebyshev2.hpp tsgRuleFejer.hpp \
		tsgRuleGaussGegenbauer.hpp tsgRuleGaussJacobi.hpp tsgRuleGaussLaguerre.hpp tsgRuleGaussHermite.hpp \
		tsgBase1DHierarchicalRule.hpp tsgRulePieceWiseLocalZero.hpp \
		tsgBaseGrid.hpp tsgTensorRule.hpp tsgGlobalGrid.hpp tsgLocalPolynomialGrid.hpp tsgFullTensorGrid.hpp \
		tsgHardcodedConstants.hpp \
		TasmanianSparseGrid.hpp \
		tasgridTestFunctions.hpp tasgridExternalTester.hpp tasgridWrapper.hpp\
		tsgRuleWavelet.hpp tsgWaveletGrid.hpp tsgSparseMatrices.hpp

LIBOBJ = tsgIndexSet.obj tsgHelperFunctions.obj \
		tsgBase1DRule.obj tsgRuleClenshawCurtis.obj tsgRuleChebyshev.obj tsgRuleGaussLegendre.obj tsgRulePieceWiseLocal.obj \
		tsgRuleChebyshevNestedTwoPoint.obj tsgRuleGaussChebyshev1.obj tsgRuleGaussChebyshev2.obj tsgRuleFejer.obj \
		tsgRuleGaussGegenbauer.obj tsgRuleGaussJacobi.obj tsgRuleGaussLaguerre.obj tsgRuleGaussHermite.obj \
		tsgBase1DHierarchicalRule.obj tsgRulePieceWiseLocalZero.obj \
		tsgBaseGrid.obj tsgTensorRule.obj tsgGlobalGrid.obj tsgLocalPolynomialGrid.obj tsgFullTensorGrid.obj \
		TasmanianSparseGrid.obj \
		tsgRuleWavelet.obj tsgWaveletGrid.obj tsgSparseMatrices.obj \
		gamma.obj
		
WROBJ = tasgrid_main.obj tasgridTestFunctions.obj tasgridExternalTester.obj tasgridWrapper.obj

OBJFILES = $(foreach file, $(LIBOBJ), $(BUILD_DIR)/$(file))
WROBJFILES = $(foreach file, $(WROBJ), $(BUILD_DIR)/$(file))

$(BUILD_DIR):
	mkdir $(BUILD_DIR)

ifeq ($(COMPILER), gcc)
LIBNAME = libtasmaniansparsegrid.so
PYLIBNAME = _py_tsg.so

$(LIBNAME): $(OBJFILES)
	$(LINK) -o $@ $(OBJFILES) $(LINKFLAGS)

$(PYLIBNAME): $(BUILD_DIR)/tsg_python.obj
	$(LINK) -o $@ $< $(LINKFLAGS) $(LIB_DIRS) $(LIBS) -L. -ltasmaniansparsegrid

$(BUILD_DIR)/tasgrid:  libtasmaniansparsegrid$(DL_SUFFIX) $(WROBJFILES)	
	$(LINK) $(LINKFLAGS_EXE) $(WROBJFILES) -o $@ -L. -ltasmaniansparsegrid

$(BUILD_DIR)/example:  libtasmaniansparsegrid$(DL_SUFFIX) $(BUILD_DIR)/example.obj	
	$(LINK) $(LINKFLAGS_EXE) -L. -ltasmaniansparsegrid $(BUILD_DIR)/example.obj -o $@
else
LIBNAME = tasmaniansparsegrid.dll
PYLIBNAME = _py_tsg.pyd

$(LIBNAME): $(OBJFILES)
	$(LINK) $(LINKFLAGS) $(OBJFILES) /OUT:$@ \
		/IMPLIB:tasmaniansparsegrid.lib \
		/MANIFESTFILE:tasmaniansparsegrid.dll.manifest

$(PYLIBNAME): $(BUILD_DIR)/tsg_python.obj
	$(LINK) $(LINKFLAGS) $(LIB_DIRS) $(LIBS) tasmaniansparsegrid.lib $< /OUT:$@
		
$(BUILD_DIR)/tasgrid.exe:  tasmaniansparsegrid$(DL_SUFFIX) $(WROBJFILES)	
	$(LINK) tasmaniansparsegrid.lib $(WROBJFILES) /OUT:$@

$(BUILD_DIR)/example.exe:  tasmaniansparsegrid$(DL_SUFFIX) $(BUILD_DIR)/example.obj	
	$(LINK) tasmaniansparsegrid.lib $(BUILD_DIR)/example.obj /OUT:$@

endif

all: $(BUILD_DIR) $(BUILD_DIR)/tasgrid$(EXE_PREFIX) $(BUILD_DIR)/example$(EXE_PREFIX) \
	$(LIBNAME) $(PYLIBNAME)
		
clean:
ifeq ($(COMPILER), gcc)
	rm -rf $(BUILD_DIR) tasmaniansparsegrid$(DL_SUFFIX)
	mkdir $(BUILD_DIR)
else
	del /q tasmaniansparsegrid$(DL_SUFFIX) tasmaniansparsegrid.lib
	rmdir /q /s $(BUILD_DIR)
	mkdir $(BUILD_DIR)
endif

