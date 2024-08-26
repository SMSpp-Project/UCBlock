##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of nc4generator                                                 #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# module name
NAME = nc4generator

# basic directory
DIR = .

# common flags
COMMON_SW = -std=c++17

# debug switches
SW_DEBUG = -g -glldb -fno-inline $(COMMON_SW)
# production switches
SW_RELEASE = -O3 $(COMMON_SW) -DNDEBUG

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

BUILD_MODE = release

default: $(BUILD_MODE)

# debug target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug: SW = $(SW_DEBUG)
debug: BUILD_MODE = debug
debug: $(DIR)/$(NAME)

# release target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

release: SW = $(SW_RELEASE)
release: BUILD_MODE = release
release: $(DIR)/$(NAME)

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(DIR)/$(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libraries + -L< libdirs >
# *H   is the list of all include files
# *INC is the -I< include directories >

# define input macros for SMS++ library makefile, then include it
# this is perhaps a bit overkill since we only need SMSTypedefs.h
SMS++SDR = ../../SMS++
include $(SMS++SDR)/lib/makefile-inc

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  # $(SMS++OBJ) no, in fact we don't need the objects

# libraries
MLIB =  $(SMS++LIB)  # ... but we do need the libs

$(DIR)/$(NAME): $(MOBJ) $(DIR)/$(NAME).o
	$(CC) -o $(DIR)/$(NAME) $(DIR)/$(NAME).o $(MLIB) $(SW)

# dependencies: every .o from its .cpp + every recursively included .h- - - -

# include directives
MINC =  $(SMS++INC)

# includes
MH =    $(SMS++H)

# compile command

$(DIR)/$(NAME).o: $(DIR)/$(NAME).cpp $(MH)
	$(CC) -c $*.cpp -o $@ $(MINC) $(SW)

# install/uninstall - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Check if the user is in the sudoers file
ifeq ($(shell id -u),0)
    prefix = /opt/smspp-project
else
    prefix = $(pwd)/smspp-project
endif
bindir = $(prefix)/$(BUILD_MODE)/bin

install: $(DIR)/$(NAME)
	install -s -D $(DIR)/$(NAME) $(DESTDIR)$(bindir)/$(NAME)

uninstall:
	rm -f $(DESTDIR)$(bindir)/$(NAME)

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: debug release install uninstall clean distclean

############################ End of makefile #################################
