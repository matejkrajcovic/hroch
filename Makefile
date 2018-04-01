# executable name
TARGET=hroch.bin

# directories
SRCDIR=src
OBJDIR=obj
BINDIR=./

SRCDIR_LIKELIHOOD=lib_likelihood
OBJDIR_LIKELIHOOD=obj_likelihood
TARGET_LIKELIHOOD=$(OBJDIR_LIKELIHOOD)/liblikelihood.a

# sources
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
INCLUDES= $(wildcard $(SRCDIR)/*.h)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

SOURCES_LIKELIHOOD = $(wildcard $(SRCDIR_LIKELIHOOD)/*.cpp)
INCLUDES_LIKELIHOOD= $(wildcard $(SRCDIR_LIKELIHOOD)/*.h)
OBJECTS_LIKELIHOOD = $(SOURCES_LIKELIHOOD:$(SRCDIR_LIKELIHOOD)/%.cpp=$(OBJDIR_LIKELIHOOD)/%.o)

# compiler and flags
CXX=g++
CXXFLAGS=-std=c++17 -static
RM=rm -f
WFLAGS=-Wall -Wextra -Wno-unused-result

# commands
all: release

build: obj $(BINDIR)/$(TARGET) outputs stats

build_likelihood: obj_likelihood $(TARGET_LIKELIHOOD)

release: CXXFLAGS+=-O2
release: build_likelihood build

debug: CXXFLAGS+=-O0 -g
debug: build_likelihood build

$(BINDIR)/$(TARGET): $(OBJECTS) $(TARGET_LIKELIHOOD)
	$(CXX) -o $@ $(OBJECTS) $(TARGET_LIKELIHOOD)
	@echo "Done."

$(TARGET_LIKELIHOOD): $(OBJECTS_LIKELIHOOD)
	$(AR) rvs $@ $(OBJECTS_LIKELIHOOD)
	@echo "Done (lib_likelihood)."

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp $(INCLUDES)
	$(CXX) $(CXXFLAGS) $(WFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully."

$(OBJECTS_LIKELIHOOD): $(OBJDIR_LIKELIHOOD)/%.o : $(SRCDIR_LIKELIHOOD)/%.cpp $(INCLUDES_LIKELIHOOD)
	$(CXX) $(CXXFLAGS) $(WFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully."

server: all
	scp -r Makefile src bio:dups/

outputs:
	mkdir -p outputs
obj:
	mkdir -p $(OBJDIR)
obj_likelihood:
	mkdir -p $(OBJDIR_LIKELIHOOD)
stats:
	mkdir -p stats

clean:
	$(RM) $(OBJECTS)
	$(RM) $(OBJECTS_LIKELIHOOD)
	$(RM) $(BINDIR)/$(TARGET)
	@echo "Clean complete."
