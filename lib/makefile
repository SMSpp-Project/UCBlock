##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of the UCBlock project                                          #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by the UCBlock project, i.e., for now the core SMS++.   #
#   These are *not* copied into $(UCBckINC): adding those -I directives to   #
#   the compile commands will have to done by whatever "main" makefile is    #
#   using this. Analogously, any external library and the corresponding      #
#   -L< libdirs > will have to be added to the final linking command by      #
#   whatever "main" makefile is using this.                                  #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by the UCBlock project.                                             #
#                                                                            #
#   Input:  $(CC)       = compiler command                                   #
#           $(SW)       = compiler options                                   #
#           $(SMS++INC) = the -I$( core SMS++ directory )                    #
#           $(SMS++OBJ) = the core SMS++ library                             #
#           $(UCBckDIR) = the directory where the source is                  #
#                                                                            #
#   Output: $(UCBckOBJ) = the final object(s) / library                      #
#           $(UCBckH)   = the .h files to include                            #
#           $(UCBckINC) = the -I$( source directory )                        #
#                                                                            #
#                                VERSION 2.00                                #
#                               12 - 02 - 2019                               #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Operations Research Group                         #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

UCBckOBJ = $(UCBckDIR)obj/UCBlock.o $(UCBckDIR)obj/NetWorkBlock.o \
	$(UCBckDIR)obj/BusNetworkBlock.o $(UCBckDIR)obj/UnitBlock.o \
	$(UCBckDIR)obj/AcadThermalUnitBlock.o \
	$(UCBckDIR)obj/AcadThemalUnitWRSolver.o \
	$(UCBckDIR)obj/AcadThemalUnitGraphSolver.o \
	$(UCBckDIR)obj/AcadThermalUnitMIPBlock.o \
	$(UCBckDIR)obj/EDPSolver.o $(UCBckDIR)obj/AcadThermalUnitMPBlock.o

UCBckINC = -I$(UCBckDIR)include/

UCBckH = $(UCBckDIR)include/UCBcklock.h $(UCBckDIR)include/NetWorkBlock.h \
	$(UCBckDIR)include/BusNetworkBlock.h $(UCBckDIR)include/UnitBlock.h \
	$(UCBckDIR)include/AcadThermalUnitBlock.h \
	$(UCBckDIR)include/AcadThemalUnitWRSolver.h \
	$(UCBckDIR)include/AcadThemalUnitGraphSolver.h \
	$(UCBckDIR)include/AcadThermalUnitMIPBlock.h \
	$(UCBckDIR)include/EDPSolver.h \
	$(UCBckDIR)include/AcadThermalUnitMPBlock.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(UCBckDIR)*.o $(UCBckDIR)*~

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(UCBckDIR)obj/UCBlock.o: $(UCBckDIR)src/UCBlock.cpp \
	$(UCBckDIR)include/UCBlock.h $(UCBckDIR)include/UnitBlock.h \
	$(UCBckDIR)include/NetWorkBlock.h $(SMS++OBJ)
	$(CC) -c $(UCBckDIR)src/UCBlock.cpp -o $@ $(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)obj/NetWorkBlock.o: $(UCBckDIR)src/NetWorkBlock.cpp \
	$(UCBckDIR)include/NetWorkBlock.h $(UCBckDIR)include/UCBlock.h \
	$(SMS++OBJ)
	$(CC) -c $(UCBckDIR)src/NetWorkBlock.cpp -o $@ $(SMS++INC) \
	$(UCBckINC) $(SW)

$(UCBckDIR)obj/BusNetworkBlock.o: $(UCBckDIR)src/BusNetworkBlock.cpp \
	$(UCBckDIR)include/BusNetworkBlock.h $(UCBckDIR)include/UCBlock.h \
	$(UCBckDIR)include/NetWorkBlock.h $(SMS++OBJ)
	$(CC) -c $(UCBckDIR)src/BusNetworkBlock.cpp -o $@ $(SMS++INC) \
	$(UCBckINC) $(SW)

$(UCBckDIR)obj/UnitBlock.o: $(UCBckDIR)src/UnitBlock.cpp \
	$(UCBckDIR)include/UnitBlock.h $(SMS++OBJ) 
	$(CC) -c $(UCBckDIR)src/UnitBlock.cpp -o $@ $(SMS++INC) \
	$(UCBckINC) $(SW)

$(UCBckDIR)obj/AcadThermalUnitBlock.o: \
	$(UCBckDIR)src/AcadThermalUnitBlock.cpp \
        $(UCBckDIR)include/AcadThermalUnitBlock.h \
	$(UCBckDIR)include/UnitBlock.h $(SMS++OBJ)
	$(CC) -c $(UCBckDIR)src/AcadThermalUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)obj/AcadThemalUnitWRSolver.o: \
	$(UCBckDIR)src/AcadThemalUnitWRSolver.cpp \
        $(UCBckDIR)include/AcadThemalUnitWRSolver.h \
	$(UCBckDIR)include/AcadThermalUnitBlock.h \
	$(UCBckDIR)include/EDPSolver.h $(SMS++OBJ)
	$(CC) -c $(UCBckDIR)src/AcadThemalUnitWRSolver.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)obj/AcadThemalUnitGraphSolver.o: \
	$(UCBckDIR)src/AcadThemalUnitGraphSolver.cpp \
        $(UCBckDIR)include/AcadThemalUnitGraphSolver.h \
	$(UCBckDIR)include/AcadThermalUnitBlock.h \
	$(UCBckDIR)include/EDPSolver.h $(SMS++OBJ)
	$(CC) -c $(UCBckDIR)src/AcadThemalUnitGraphSolver.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)obj/AcadThermalUnitMIPBlock.o: \
	$(UCBckDIR)src/AcadThermalUnitMIPBlock.cpp \
        $(UCBckDIR)include/AcadThermalUnitMIPBlock.h \
	$(UCBckDIR)include/UnitBlock.h $(SMS++OBJ)
	$(CC) -c $(UCBckDIR)src/AcadThermalUnitMIPBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)obj/EDPSolver.o: $(UCBckDIR)src/EDPSolver.cpp \
        $(UCBckDIR)include/EDPSolver.h $(UCBckDIR)AcadThermalUnitBlock.h \
	$(UCBckDIR)include/UnitBlock.h $(SMS++OBJ)
	$(CC) -c $(UCBckDIR)src/EDPSolver.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)obj/AcadThermalUnitMPBlock.o: \
	$(UCBckDIR)src/AcadThermalUnitMPBlock.cpp \
        $(UCBckDIR)include/AcadThermalUnitMIPBlock.h \
	$(UCBckDIR)include/UnitBlock.h $(SMS++OBJ)
	$(CC) -c $(UCBckDIR)src/AcadThermalUnitMPBlock.cpp -o \
	$@ $(SMS++INC) $(UCBckINC) $(SW)

########################## End of makefile ###################################
