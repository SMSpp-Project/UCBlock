##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test (for core SMS++)                                        #
#                                                                            #
#                                VERSION 2.00                                #
#                               20 - 04 - 2005                               #
#                                                                            #
#                              Implementation by:                            #
#                                                                            #
#                              Antonio Frangioni                             #
#                           Operations Research Group                        #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# module name
NAME = test

# basic directory
DIR = ./

# debug switches
SW =-g  \
     -fmax-errors=1 -std=c++14


#SW =  -g -Wall -Werror -Wno-unused-local-typedefs -Wno-sign-compare \
 #    -fmax-errors=1 -std=c++17 -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD \
  #   -D_GLIBCXX_USE_CXX11_ABI=0 -lstdc++fs



# production switches
#SW = -O3 -std=c++14

# compiler
CC = g++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: $(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)*.o $(DIR)*~ $(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libreries + -L<libdirs>
# *H   is the list of all include files
# *INC is the -I<include directories>

# core SMS++
SMSDIR = ../../SMS++/
SMSINC = -I$(SMSDIR)/include/
SMSH = $(SMSDIR)/include/*.h
#SMSLIB = $(SMSDIR)/lib/
SMSOBJ = $(SMSDIR)/obj/*.o
#include $(SMSDIR)/test/makefile

# UCBlock
UCBDIR = ../UCBLock/
include $(UCBDIR)makefile

# CPLEX - - - - - - - - - - - - - - - - - - - - - - - -
#CPLXDIR       = /usr/cplex128/cplex/examples/src/c/
#CPLEXDIR      = /usr/cplex128/cplex/include
#CONCERTDIR    = /usr/cplex128/concert/include
#CPLEXLIBDIR   = /usr/cplex128/cplex/lib/x86-64_osx/static_pic
#CONCERTLIBDIR = /usr/cplex128/concert/lib/x86-64_osx/static_pic

CPLXDIR       = /Applications/CPLEX_Studio128/cplex/examples/src/c/
CPLEXDIR      = /Applications/CPLEX_Studio128/cplex/include
CONCERTDIR    = /Applications/CPLEX_Studio128/concert/include
CPLEXLIBDIR   = /Application/CPLEX_Studio128/cplex/lib/x86-64_osx/static_pic
CONCERTLIBDIR = /Applications/CPLEX_Studio128/concert/lib/x86-64_osx/static_pic


CPLEXLIB = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) -lilocplex -lconcert -lcplex -lm -lpthread

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -


# object files
MOBJ =  $(SMSOBJ) $(UCBOBJ)

# libreries
MLIB =  $(SMSLIB) $(CPLEXLIB)



#$(NAME): $(MOBJ) $(DIR)test.o $(CPLXDIR)check.o
#	$(CC) -o $(NAME) $(DIR)test.o $(MOBJ) $(CPLXDIR)check.o $(MLIB) $(SW)

$(NAME): $(MOBJ) $(DIR)test.o 
	$(CC) -o $(NAME) $(DIR)test.o $(MOBJ) $(MLIB) $(SW)

# dependencies: every .o from its .C + every recursively included .h- - - - -

# include dirs

CPLEXINC = -I$(CPLEXDIR)  -I$(CONCERTDIR) 

MINC =  $(SMSINC) $(UCBINC) $(CPLEXINC)

# includes
MH =    $(SMSH) $(UCBH)



$(DIR)test.o: $(DIR)test.cpp $(MH)
	$(CC) -c $*.cpp -o $@ $(MINC)  $(SW)

#$(CPLXDIR)check.o: $(CPLXDIR)check.c 
#	$(CC) -c $*.c -o $@ $(CPLEXINC)  $(SW)


############################ End of makefile #################################
