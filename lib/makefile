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
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

UCBckOBJ = $(UCBckDIR)/obj/BatteryUnitBlock.o \
	$(UCBckDIR)/obj/DCNetworkBlock.o \
	$(UCBckDIR)/obj/ECNetworkBlock.o \
	$(UCBckDIR)/obj/HydroSystemUnitBlock.o \
	$(UCBckDIR)/obj/HydroUnitBlock.o \
	$(UCBckDIR)/obj/IntermittentUnitBlock.o \
	$(UCBckDIR)/obj/NetworkBlock.o \
	$(UCBckDIR)/obj/NuclearUnitBlock.o \
	$(UCBckDIR)/obj/SlackUnitBlock.o \
	$(UCBckDIR)/obj/ThermalUnitBlock.o \
	$(UCBckDIR)/obj/ThermalUnitDPSolver.o \
	$(UCBckDIR)/obj/UCBlock.o \
	$(UCBckDIR)/obj/UnitBlock.o
#	$(UCBckDIR)/obj/HeatBlock.o

UCBckINC = -I$(UCBckDIR)/include

UCBckH = $(UCBckDIR)/include/BatteryUnitBlock.h \
	$(UCBckDIR)/include/DCNetworkBlock.h \
	$(UCBckDIR)/include/ECNetworkBlock.h \
	$(UCBckDIR)/include/HydroSystemUnitBlock.h \
	$(UCBckDIR)/include/HydroUnitBlock.h \
	$(UCBckDIR)/include/IntermittentUnitBlock.h \
	$(UCBckDIR)/include/NetworkBlock.h \
	$(UCBckDIR)/include/NuclearUnitBlock.h \
	$(UCBckDIR)/include/SlackUnitBlock.h \
	$(UCBckDIR)/include/ThermalUnitBlock.h \
	$(UCBckDIR)/include/ThermalUnitDPSolver.h \
	$(UCBckDIR)/include/UCBlock.h \
	$(UCBckDIR)/include/UnitBlock.h
#	$(UCBckDIR)/include/HeatBlock.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(UCBckDIR)/*.o $(UCBckDIR)/*~

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(UCBckDIR)/obj/BatteryUnitBlock.o: $(UCBckDIR)/src/BatteryUnitBlock.cpp \
        $(UCBckDIR)/include/BatteryUnitBlock.h $(UCBckDIR)/include/UnitBlock.h \
	$(SMS++OBJ) 
	$(CC) -c $(UCBckDIR)/src/BatteryUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)/obj/DCNetworkBlock.o: $(UCBckDIR)/src/DCNetworkBlock.cpp \
	$(UCBckDIR)/include/DCNetworkBlock.h $(UCBckDIR)/include/UCBlock.h \
	$(UCBckDIR)/include/NetworkBlock.h $(SMS++OBJ)
	$(CC) -c $(UCBckDIR)/src/DCNetworkBlock.cpp -o $@ $(SMS++INC) \
	$(UCBckINC) $(SW)

$(UCBckDIR)/obj/ECNetworkBlock.o: $(UCBckDIR)/src/ECNetworkBlock.cpp \
	$(UCBckDIR)/include/ECNetworkBlock.h $(UCBckDIR)/include/UCBlock.h \
	$(UCBckDIR)/include/NetworkBlock.h $(SMS++OBJ)
	$(CC) -c $(UCBckDIR)/src/ECNetworkBlock.cpp -o $@ $(SMS++INC) \
	$(UCBckINC) $(SW)

#$(UCBckDIR)/obj/HeatBlock.o: $(UCBckDIR)/src/HeatBlock.cpp \
#	$(UCBckDIR)/include/HeatBlock.h $(UCBckDIR)/include/UnitBlock.h \
#	$(SMS++OBJ) 
#	$(CC) -c $(UCBckDIR)/src/HeatBlock.cpp -o $@ \
#	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)/obj/HydroSystemUnitBlock.o: $(UCBckDIR)/src/HydroSystemUnitBlock.cpp \
        $(UCBckDIR)/include/HydroSystemUnitBlock.h \
	$(UCBckDIR)/include/HydroUnitBlock.h $(UCBckDIR)/include/UnitBlock.h \
	$(SMS++OBJ) 
	$(CC) -c $(UCBckDIR)/src/HydroSystemUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)/obj/HydroUnitBlock.o: $(UCBckDIR)/src/HydroUnitBlock.cpp \
	$(UCBckDIR)/include/HydroUnitBlock.h $(UCBckDIR)/include/UnitBlock.h \
	$(SMS++OBJ) 
	$(CC) -c $(UCBckDIR)/src/HydroUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)/obj/IntermittentUnitBlock.o: \
	$(UCBckDIR)/src/IntermittentUnitBlock.cpp \
	$(UCBckDIR)/include/IntermittentUnitBlock.h \
	$(UCBckDIR)/include/UnitBlock.h $(SMS++OBJ) 
	$(CC) -c $(UCBckDIR)/src/IntermittentUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)/obj/NetworkBlock.o: $(UCBckDIR)/src/NetworkBlock.cpp \
	$(UCBckDIR)/include/NetworkBlock.h $(UCBckDIR)/include/UCBlock.h \
	$(SMS++OBJ)
	$(CC) -c $(UCBckDIR)/src/NetworkBlock.cpp -o $@ $(SMS++INC) \
	$(UCBckINC) $(SW)

$(UCBckDIR)/obj/NuclearUnitBlock.o: $(UCBckDIR)/src/NuclearUnitBlock.cpp \
        $(UCBckDIR)/include/NuclearUnitBlock.h \
	$(UCBckDIR)/include/ThermalUnitBlock.h \
	$(UCBckDIR)/include/UnitBlock.h $(SMS++OBJ) 
	$(CC) -c $(UCBckDIR)/src/NuclearUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)/obj/SlackUnitBlock.o: $(UCBckDIR)/src/SlackUnitBlock.cpp \
	$(UCBckDIR)/include/SlackUnitBlock.h $(UCBckDIR)/include/UnitBlock.h \
	$(SMS++OBJ) 
	$(CC) -c $(UCBckDIR)/src/SlackUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)/obj/ThermalUnitBlock.o: $(UCBckDIR)/src/ThermalUnitBlock.cpp \
        $(UCBckDIR)/include/ThermalUnitBlock.h $(UCBckDIR)/include/UnitBlock.h \
	$(SMS++OBJ) 
	$(CC) -c $(UCBckDIR)/src/ThermalUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)/obj/ThermalUnitDPSolver.o: $(UCBckDIR)/src/ThermalUnitDPSolver.cpp \
        $(UCBckDIR)/include/ThermalUnitDPSolver.h $(SMS++OBJ)
	$(CC) -c $(UCBckDIR)/src/ThermalUnitDPSolver.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)/obj/UCBlock.o: $(UCBckDIR)/src/UCBlock.cpp \
	$(UCBckDIR)/include/UCBlock.h $(UCBckDIR)/include/UnitBlock.h \
	$(UCBckDIR)/include/NetworkBlock.h $(SMS++OBJ)
	$(CC) -c $(UCBckDIR)/src/UCBlock.cpp -o $@ $(SMS++INC) $(UCBckINC) $(SW)

$(UCBckDIR)/obj/UnitBlock.o: $(UCBckDIR)/src/UnitBlock.cpp \
	$(UCBckDIR)/include/UnitBlock.h $(SMS++OBJ) 
	$(CC) -c $(UCBckDIR)/src/UnitBlock.cpp -o $@ $(SMS++INC) $(UCBckINC) $(SW)

########################## End of makefile ###################################
