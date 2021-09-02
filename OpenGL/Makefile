####
#### determine what platform: linux, osx or windows
####
ifeq ($(OS),Windows_NT)
  detected_OS := Windows
else
  detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')
endif

ifeq ($(detected_OS),Windows)
    PLATFORM = WINDOWS
endif
ifeq ($(detected_OS),Darwin)
    PLATFORM = OSX
endif
ifeq ($(detected_OS),Linux)
    PLATFORM = LINUX
endif



# create directory command
MKDIR   = mkdir -p


# Set build directory
BUILDDIR = build


# Set source directory and source files extension
SRCDIR = src
SRCEXT = cpp


# Set executable name and directory
TARGETDIR = bin
TARGET = main


# set library directory
LIBDIR = libs


# create compiler flags
CXX       = g++ -std=c++11 -Wall --pedantic
INCFLAGS  = -I ./$(LIBDIR)/build/include/
LIBFLAGS  = -L ./$(LIBDIR)/build/lib/
LIBS      = -lglad -lglfw3

# # AntTweakBar
# INCFLAGS  += -I ./$(LIBDIR)/AntTweakBar/include/
# LIBFLAGS  += -L ./$(LIBDIR)/AntTweakBar/lib/
# LIBS      += -l AntTweakBar

# platform dependent part
ifeq ($(PLATFORM),LINUX)
	# linux
	LIBS += -lEGL -lXext -ldl -pthread -lGL -lX11 -lXi -lXrandr -lXxf86vm -lXinerama -lXcursor -lrt -lm
endif
ifeq ($(PLATFORM),OSX)
	# mac
	LIBS += -framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo
endif


####
#### ------------- HERE FOLLOW TARGETS -------------
####

# Explore source tree
SRCS := $(shell find $(SRCDIR) -type f -name '*.$(SRCEXT)')


# Create build tree structure: src/../*.cpp -> build/../*.o
OBJS := $(subst $(SRCDIR),$(BUILDDIR),$(SRCS:.$(SRCEXT)=.o))


# Define dependencies files for all objects
DEPS = $(OBJS:.o=.d)


# Indicate to make which targets are not files
.PHONY: all clean setup platform


# The default target triggered by the "make" command
# First the setup rule is run, then the build rule
all: setup $(TARGETDIR)/$(TARGET)


# The target creating all necessary directories for the build process
# Triggered by "make setup"
setup:
	@$(MKDIR) $(dir $(OBJS))
	@$(MKDIR) $(TARGETDIR)


# The build process target, which is a file
$(TARGETDIR)/$(TARGET): $(OBJS)
	@echo linking: $@
	@$(CXX) $(LIBFLAGS) -o $@ $(OBJS) $(LIBS)


# Defines the function which will generate the rule for each object (.o file)
define generateRules =
$(1): $(subst .o,.cpp,$(subst $(BUILDDIR),$(SRCDIR),$(1)))
	@echo compiling: $$<
	@$$(CXX) $(INCFLAGS) -c -o $$@ $$< -MMD
endef

# Generates the rules for each .o file
$(foreach obj,$(OBJS),$(eval $(call generateRules,$(obj))))


# Removes build tree and the binary file
# Triggered by "make clean"
clean:
	@$(RM) -rf $(TARGETDIR)/$(TARGET)
	@$(RM) -fr $(BUILDDIR)

# rebuild rule
rebuild: clean all

run: all
	@echo ------------[$(TARGETDIR)/$(TARGET)]------------
	@$(TARGETDIR)/$(TARGET)

# outputs what platform
platform:
	@echo $(PLATFORM)

# Include dependencies
-include $(DEPS)