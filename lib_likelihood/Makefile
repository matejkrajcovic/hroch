# executable name
TARGET=pramoduhis.bin

# directories
SRCDIR=source
OBJDIR=obj
BINDIR=./

# sources
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
INCLUDES= $(wildcard $(SRCDIR)/*.h)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# compiler and flags

CXX=g++
CXXFLAGS=-O2 -std=gnu++11 -static
RM=rm
WFLAGS=-Wall -Wextra -Wno-unused-result -g -pg -static

# commands

all: $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(OBJECTS)
	$(CXX) -o $@ $(WFLAGS) $(OBJECTS)
	@echo "Done."

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(WFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully."

server: all
	scp $(BINDIR)/$(TARGET) compbio.fmph.uniba.sk:Bakalarka/

clean:
	$(RM) $(OBJECTS)
	@echo "Clean complete."

remove: clean
	$(RM) $(BINDIR)/$(TARGET)
	@echo "Remove complete."
