# executable name
TARGET=hroch.bin

# directories
SRCDIR=src
OBJDIR=obj
BINDIR=./

# sources
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
INCLUDES= $(wildcard $(SRCDIR)/*.h)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# compiler and flags

CXX=g++
CXXFLAGS=-std=gnu++11 -O2 -static
RM=rm
WFLAGS=-Wall -Wextra -Wno-unused-result
#-g -pg -static

# commands

all: obj $(BINDIR)/$(TARGET) outputs stats

$(BINDIR)/$(TARGET): $(OBJECTS)
	$(CXX) -o $@ $(WFLAGS) $(OBJECTS)
	@echo "Done."

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp $(INCLUDES)
	$(CXX) $(CXXFLAGS) $(WFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully."

server: all
	scp -r Makefile src bio:dups/

outputs:
	mkdir -p outputs
obj:
	mkdir -p obj
stats:
	mkdir -p stats

clean:
	$(RM) $(OBJECTS)
	@echo "Clean complete."

remove: clean
	$(RM) $(BINDIR)/$(TARGET)
	@echo "Remove complete."
