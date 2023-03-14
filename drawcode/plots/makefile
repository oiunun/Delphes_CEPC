CC         =  gcc
CPP        =  g++
CXX        =  g++
LINK       =  g++
LD         =  g++
# 
CCFLAGS  = -w -O2 -Wall 
DEBUG    = -g    
PREPRO   =
LFLAGS   =   $(shell root-config --cflags)
INCPATH  = -I$(shell root-config --incdir)  -I/usr/X11R6/include -I./
LIBS     =   $(shell root-config --libs)
GLIBS    =   $(shell root-config --glibs)
#
AR       = ar cqs
TAR      = tar -cf
GZIP     = gzip -9f
COPY     = cp -f
COPY_FILE= $(COPY)
COPY_DIR = $(COPY) -r
INSTALL_FILE= $(COPY_FILE)
INSTALL_DIR = $(COPY_DIR)
DEL_FILE = rm -f
SYMLINK  = ln -sf
DEL_DIR  = rmdir
MOVE     = mv -f
CHK_DIR_EXISTS= test -d
MKDIR    = mkdir -p
####### Output directory
OBJECTS_DIR = ./
####### Files
HEADERS = cepcplotstyle.h     
SOURCES = plots.C  cepcplotstyle.C     
OBJECTS = plots.o  cepcplotstyle.o     
FORMS    = 
DESTDIR  = 
TARGET   = plots

####### Implicit rules
# creating extra dictionary sources
%Cint.C: %.h LinkDef.h
	@echo "Generating dictionary ..."
	@rootcint -f %Cint.C -c %.h %Cint.h
# compiling c++ code, for root
%.o: %.C
	@echo $(@) $^
	$(CPP) -c $(CCFLAGS) $(INCPATH) $(DEFINES) $(DEBUG)\
	-o $@ $<
first: all
####### Build rules
all: $(TARGET) 
$(TARGET):  $(OBJECTS) 
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS)  $(OBJCOMP) $(LIBS)
clean:
	-$(DEL_FILE) $(OBJECTS) $(OBJECTS1)
	-$(DEL_FILE) *~ core *.core
	-$(DEL_FILE) *Cint.?
####### Sub-libraries
distclean: clean
	-$(DEL_FILE) $(TARGET) $(TARGET)
FORCE:
####### Compile
 plots.o           : plots.C      cepcplotstyle.h     
 cepcplotsyle.o    : cepcplotstyle.C   cepcplotstyle.h     
####### Install
install:  
uninstall:  
