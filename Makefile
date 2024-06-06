PACS_ROOT = /home/irene/PACS/pacs-examples/Examples
CXX = mpic++
CXXFLAGS = -std=c++20
CPPFLAGS=-fopenmp -O3 -Wall -pedantic -I$(PACS_ROOT)/include 
SRC_DIR = src
INCLUDE_DIR = include 
TARGET = main
DOXYFILE=Doxyfile
LDFLAGS = -L$(PACS_ROOT)/lib
LIBS  = -lmuparser -lgomp
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,%.o,$(SRCS))
INCLUDES:= -I$(INCLUDE_DIR)

.PHONY: all clean
all: $(TARGET)
%.o: $(SRC_DIR)/%.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $<


$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< $(LIBS) -o $@

clean:
	$(RM) $(OBJS) $(TARGET) output/*.vtk

distclean: clean
	$(RM) *~

doc:
	doxygen $(DOXYFILE)
