##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test (for the UCBlock project)                               #
#                                                                            #
#                                VERSION 1.00                                #
#                               12 - 02 - 2019                               #
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
SW = -g -std=c++17

# production switches
#SW = -O3 -DNDEBUG -std=c++14

# compiler
#CC = g++-8
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

all default: $(NAME)

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

# UCBlock project
UCBckDIR = ../
include $(UCBckDIR)lib/makefile-inc

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(UCBckOBJ)

# libreries
MLIB =  $(UCBckLIB)

$(NAME): $(MOBJ) $(DIR)$(NAME).o
	$(CC) -o $(NAME) $(DIR)$(NAME).o $(MOBJ) $(MLIB) $(SW)

# dependencies: every .o from its .C + every recursively included .h- - - - -

# include directives
MINC =  $(UCBckINC)

# includes
MH =    $(UCBckH)

# compile command

$(DIR)$(NAME).o: $(DIR)$(NAME).cpp $(MH)
	$(CC) -c $*.cpp -o $@ $(MINC) $(SW)

############################ End of makefile ################################

