# Compiler and flags
CXX = cl

CXXFLAGS = /std:c++17 /Zi /O2 /W3 /EHsc /MD /DEBUG:FULL /JMC
INCLUDES = /I"include" /I"D:/gurobi1202/win64/include"
LIBS = /link /LIBPATH:"D:/gurobi1202/win64/lib" gurobi_c++md2017.lib gurobi120.lib

# Directories
SRCDIR = src
BUILDDIR = build
TARGET = $(BUILDDIR)\bmf_solver.exe

# Source and object files
SRC = $(wildcard $(SRCDIR)/*.cpp)
OBJ = $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.obj,$(SRC))

# Default target
all: $(TARGET)

# Build target
$(TARGET): main.cpp $(OBJ)
	@if not exist $(BUILDDIR) mkdir $(BUILDDIR)
	$(CXX) main.cpp $(OBJ) $(CXXFLAGS) $(INCLUDES) $(LIBS) /OUT:$(TARGET)

# Compile source files to object files
$(BUILDDIR)/%.obj: $(SRCDIR)/%.cpp
	@if not exist $(BUILDDIR) mkdir $(BUILDDIR)
	$(CXX) /c $(CXXFLAGS) $(INCLUDES) $< /Fo$@

# Clean build
clean:
	@if exist $(BUILDDIR) del /Q $(BUILDDIR)\*.obj
	@if exist $(TARGET) del /Q $(TARGET)

.PHONY: all clean
