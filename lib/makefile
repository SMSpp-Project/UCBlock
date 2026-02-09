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

# internal macros- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

UCBckINCdir = $(UCBckDIR)/include
UCBckSRCdir = $(UCBckDIR)/src
UCBckOdir = $(UCBckDIR)/obj

# macros to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

UCBckOBJ = $(UCBckOdir)/ACNetworkBlock.o \
	$(UCBckOdir)/BatteryUnitBlock.o \
	$(UCBckOdir)/DesignNetworkBlock.o \
	$(UCBckOdir)/DCNetworkBlock.o \
	$(UCBckOdir)/ECNetworkBlock.o \
	$(UCBckOdir)/HydroSystemUnitBlock.o \
	$(UCBckOdir)/HydroUnitBlock.o \
	$(UCBckOdir)/IntermittentUnitBlock.o \
	$(UCBckOdir)/NetworkBlock.o \
	$(UCBckOdir)/NuclearUnitBlock.o \
	$(UCBckOdir)/SlackUnitBlock.o \
	$(UCBckOdir)/ThermalUnitBlock.o \
	$(UCBckOdir)/ThermalUnitDPSolver.o \
	$(UCBckOdir)/UCBlock.o \
	$(UCBckOdir)/UnitBlock.o

UCBckINC = -I$(UCBckINCdir)/

UCBckH = $(UCBckINCdir)/BatteryUnitBlock.h \
	$(UCBckINCdir)/DesignNetworkBlock.h \
	$(UCBckINCdir)/DCNetworkBlock.h \
	$(UCBckINCdir)/ECNetworkBlock.h \
	$(UCBckINCdir)/HydroSystemUnitBlock.h \
	$(UCBckINCdir)/HydroUnitBlock.h \
	$(UCBckINCdir)/IntermittentUnitBlock.h \
	$(UCBckINCdir)/NetworkBlock.h \
	$(UCBckINCdir)/NuclearUnitBlock.h \
	$(UCBckINCdir)/SlackUnitBlock.h \
	$(UCBckINCdir)/ThermalUnitBlock.h \
	$(UCBckINCdir)/ThermalUnitDPSolver.h \
	$(UCBckINCdir)/UCBlock.h \
	$(UCBckINCdir)/UnitBlock.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(UCBckDIR)/*.o $(UCBckDIR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(UCBckOdir)/ACNetworkBlock.o: $(UCBckSRCdir)/ACNetworkBlock.cpp \
	$(UCBckINCdir)/ACNetworkBlock.h $(UCBckINCdir)/DCNetworkBlock.h \
	$(UCBckINCdir)/UCBlock.h $(UCBckINCdir)/NetworkBlock.h $(SMS++OBJ)
	$(CC) -c $(UCBckSRCdir)/ACNetworkBlock.cpp -o $@ $(SMS++INC) \
	$(UCBckINC) $(SW)

$(UCBckOdir)/BatteryUnitBlock.o: $(UCBckSRCdir)/BatteryUnitBlock.cpp \
        $(UCBckINCdir)/BatteryUnitBlock.h $(UCBckINCdir)/UnitBlock.h \
	$(SMS++OBJ) 
	$(CC) -c $(UCBckSRCdir)/BatteryUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckOdir)/DesignNetworkBlock.o: $(UCBckSRCdir)/DesignNetworkBlock.cpp \
	$(UCBckINCdir)/DesignNetworkBlock.h $(UCBckINCdir)/UCBlock.h \
	$(UCBckINCdir)/NetworkBlock.h $(UCBckINCdir)/DCNetworkBlock.h \
	$(SMS++OBJ)
	$(CC) -c $(UCBckSRCdir)/DesignNetworkBlock.cpp -o $@ $(SMS++INC) \
	$(UCBckINC) $(SW)

$(UCBckOdir)/DCNetworkBlock.o: $(UCBckSRCdir)/DCNetworkBlock.cpp \
	$(UCBckINCdir)/DCNetworkBlock.h $(UCBckINCdir)/UCBlock.h \
	$(UCBckINCdir)/NetworkBlock.h $(SMS++OBJ)
	$(CC) -c $(UCBckSRCdir)/DCNetworkBlock.cpp -o $@ $(SMS++INC) \
	$(UCBckINC) $(SW)

$(UCBckOdir)/ECNetworkBlock.o: $(UCBckSRCdir)/ECNetworkBlock.cpp \
	$(UCBckINCdir)/ECNetworkBlock.h $(UCBckINCdir)/UCBlock.h \
	$(UCBckINCdir)/NetworkBlock.h $(SMS++OBJ)
	$(CC) -c $(UCBckSRCdir)/ECNetworkBlock.cpp -o $@ $(SMS++INC) \
	$(UCBckINC) $(SW)

$(UCBckOdir)/HydroSystemUnitBlock.o: \
	$(UCBckSRCdir)/HydroSystemUnitBlock.cpp \
        $(UCBckINCdir)/HydroSystemUnitBlock.h \
	$(UCBckINCdir)/HydroUnitBlock.h $(UCBckINCdir)/UnitBlock.h \
	$(SMS++OBJ) 
	$(CC) -c $(UCBckSRCdir)/HydroSystemUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckOdir)/HydroUnitBlock.o: $(UCBckSRCdir)/HydroUnitBlock.cpp \
	$(UCBckINCdir)/HydroUnitBlock.h $(UCBckINCdir)/UnitBlock.h \
	$(SMS++OBJ) 
	$(CC) -c $(UCBckSRCdir)/HydroUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckOdir)/IntermittentUnitBlock.o: \
	$(UCBckSRCdir)/IntermittentUnitBlock.cpp \
	$(UCBckINCdir)/IntermittentUnitBlock.h $(UCBckINCdir)/UnitBlock.h \
	$(SMS++OBJ) 
	$(CC) -c $(UCBckSRCdir)/IntermittentUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckOdir)/NetworkBlock.o: $(UCBckSRCdir)/NetworkBlock.cpp \
	$(UCBckINCdir)/NetworkBlock.h $(UCBckINCdir)/UCBlock.h $(SMS++OBJ)
	$(CC) -c $(UCBckSRCdir)/NetworkBlock.cpp -o $@ $(SMS++INC) \
	$(UCBckINC) $(SW)

$(UCBckOdir)/NuclearUnitBlock.o: $(UCBckSRCdir)/NuclearUnitBlock.cpp \
        $(UCBckINCdir)/NuclearUnitBlock.h $(UCBckINCdir)/ThermalUnitBlock.h \
	$(UCBckINCdir)/UnitBlock.h $(SMS++OBJ) 
	$(CC) -c $(UCBckSRCdir)/NuclearUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckOdir)/SlackUnitBlock.o: $(UCBckSRCdir)/SlackUnitBlock.cpp \
	$(UCBckINCdir)/SlackUnitBlock.h $(UCBckINCdir)/UnitBlock.h \
	$(SMS++OBJ) 
	$(CC) -c $(UCBckSRCdir)/SlackUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckOdir)/ThermalUnitBlock.o: $(UCBckSRCdir)/ThermalUnitBlock.cpp \
        $(UCBckINCdir)/ThermalUnitBlock.h $(UCBckINCdir)/UnitBlock.h \
	$(SMS++OBJ) 
	$(CC) -c $(UCBckSRCdir)/ThermalUnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckOdir)/ThermalUnitDPSolver.o: \
	$(UCBckSRCdir)/ThermalUnitDPSolver.cpp \
        $(UCBckINCdir)/ThermalUnitDPSolver.h $(SMS++OBJ)
	$(CC) -c $(UCBckSRCdir)/ThermalUnitDPSolver.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckOdir)/UCBlock.o: $(UCBckSRCdir)/UCBlock.cpp \
	$(UCBckINCdir)/UCBlock.h $(UCBckINCdir)/UnitBlock.h \
	$(UCBckINCdir)/NetworkBlock.h $(SMS++OBJ)
	$(CC) -c $(UCBckSRCdir)/UCBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

$(UCBckOdir)/UnitBlock.o: $(UCBckSRCdir)/UnitBlock.cpp \
	$(UCBckINCdir)/UnitBlock.h $(SMS++OBJ) 
	$(CC) -c $(UCBckSRCdir)/UnitBlock.cpp -o $@ \
	$(SMS++INC) $(UCBckINC) $(SW)

########################## End of makefile ###################################
