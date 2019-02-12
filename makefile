##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of the UCBlock                                                  #
#                                                                            #
#   Input:  $(CC)     = compiler command                                     #
#           $(SW)     = compiler options                                     #
#           $(UCBDIR) = the directory where the core UCBlock source is       #
#                                                                            #
#   Output: $(UCBOBJ) = the final object(s) / library                        #
#           $(UCBLIB) = external libreries + -L<libdirs>                     #
#           $(UCBH)   = the .h files to include for core SMS++               #
#           $(UCBINC) = the -I$( core SMS++ directory)                       #
#                                                                            #
#                                VERSION 1.00                                #
#                               04 - 07 - 2016                               #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Operations Research Group                         #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(UCBDIR)*.o $(UCBDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

UCBLIB = -L/usr/local/lib -lboost_serialization
#UCBLIB = 

UCBOBJ =$(UCBDIR)UCBlock.o \
	$(UCBDIR)NetWorkBlock.o $(UCBDIR)BusNetworkBlock.o  \
	$(UCBDIR)UnitBlock.o $(UCBDIR)AcadThermalUnitBlock.o \
	$(UCBDIR)AcadThemalUnitWRSolver.o $(UCBDIR)AcadThemalUnitGraphSolver.o \
    $(UCBDIR)AcadThermalUnitMIPBlock.o $(UCBDIR)EDPSolver.o $(UCBDIR)AcadThermalUnitMPBlock.o

UCBINC = -I$(UCBDIR)

UCBH   =$(UCBDIR)UCBlock.h \
	$(UCBDIR)NetWorkBlock.h $(UCBDIR)BusNetworkBlock.h \
	$(UCBDIR)UnitBlock.h $(UCBDIR)AcadThermalUnitBlock.h \
    $(UCBDIR)AcadThemalUnitWRSolver.h $(UCBDIR)AcadThemalUnitGraphSolver.h\
    $(UCBDIR)AcadThermalUnitMIPBlock.h $(UCBDIR)EDPSolver.h $(UCBDIR)AcadThermalUnitMPBlock.h


# dependencies: every .o from its .C + every recursively included .h- - - - -

$(UCBDIR)UCBlock.o: $(UCBDIR)UCBlock.cpp $(UCBDIR)UCBlock.h \
        $(UCBDIR)UnitBlock.h $(UCBDIR)NetWorkBlock.h $(SMSDIR)Block.h
	$(CC) -c $*.cpp -o $@ $(SMSINC) $(UCBINC) $(SW)

$(UCBDIR)NetWorkBlock.o: $(UCBDIR)NetWorkBlock.cpp $(UCBDIR)NetWorkBlock.h \
        $(UCBDIR)UCBlock.h $(SMSDIR)Block.h 
	$(CC) -c $*.cpp -o $@ $(SMSINC) $(UCBINC) $(SW)

$(UCBDIR)BusNetworkBlock.o: $(UCBDIR)BusNetworkBlock.cpp $(UCBDIR)BusNetworkBlock.h  $(UCBDIR)UnitBlock.h\
        $(UCBDIR)UCBlock.h $(UCBDIR)NetWorkBlock.h $(SMSDIR)LinearConstraint.h 
	$(CC) -c $*.cpp -o $@ $(SMSINC) $(UCBINC) $(SW)

$(UCBDIR)UnitBlock.o: $(UCBDIR)UnitBlock.cpp $(UCBDIR)UnitBlock.h \
        $(UCBDIR)UCBlock.h $(SMSDIR)Block.h $(SMSDIR)ColVariable.h 
	$(CC) -c $*.cpp -o $@ $(SMSINC) $(UCBINC) $(SW)

$(UCBDIR)AcadThermalUnitBlock.o: $(UCBDIR)AcadThermalUnitBlock.cpp \
        $(UCBDIR)AcadThermalUnitBlock.h $(UCBDIR)UCBlock.h \
        $(SMSDIR)LinearConstraint.h $(SMSDIR)DQuadObjectiveFunction.h \
	$(UCBDIR)UnitBlock.h
	$(CC) -c $*.cpp -o $@ $(SMSINC) $(UCBINC) $(SW)

$(UCBDIR)AcadThemalUnitWRSolver.o: $(UCBDIR)AcadThemalUnitWRSolver.cpp \
        $(UCBDIR)AcadThemalUnitWRSolver.h $(UCBDIR)AcadThermalUnitBlock.h \
	$(SMSDIR)Block.h $(UCBDIR)UCBlock.h $(SMSDIR)Solver.h \
        $(SMSDIR)SMSTypedefs.h $(UCBDIR)EDPSolver.h 
	$(CC) -c $*.cpp -o $@ $(SMSINC) $(UCBINC) $(SW)

$(UCBDIR)AcadThemalUnitGraphSolver.o: $(UCBDIR)AcadThemalUnitGraphSolver.cpp \
        $(UCBDIR)AcadThemalUnitGraphSolver.h $(UCBDIR)AcadThermalUnitBlock.h \
	$(SMSDIR)Block.h $(UCBDIR)UCBlock.h $(SMSDIR)Solver.h \
        $(SMSDIR)SMSTypedefs.h $(UCBDIR)EDPSolver.h 
	$(CC) -c $*.cpp -o $@ $(SMSINC) $(UCBINC) $(SW)

$(UCBDIR)AcadThermalUnitMIPBlock.o: $(UCBDIR)AcadThermalUnitMIPBlock.cpp \
        $(UCBDIR)AcadThermalUnitMIPBlock.h $(SMSDIR)LinearConstraint.h \
        $(SMSDIR)DQuadObjectiveFunction.h $(UCBDIR)UCBlock.h
	$(CC) -c $*.cpp -o $@ $(SMSINC) $(UCBINC) $(SW)

$(UCBDIR)EDPSolver.o: $(UCBDIR)EDPSolver.cpp \
        $(UCBDIR)EDPSolver.h $(UCBDIR)AcadThermalUnitBlock.h \
	$(UCBDIR)UCBlock.h $(SMSDIR)SMSTypedefs.h
	$(CC) -c $*.cpp -o $@ $(SMSINC) $(UCBINC) $(SW)

$(UCBDIR)AcadThermalUnitMPBlock.o: $(UCBDIR)AcadThermalUnitMPBlock.cpp \
        $(UCBDIR)AcadThermalUnitMIPBlock.h $(SMSDIR)LinearConstraint.h \
        $(SMSDIR)DQuadObjectiveFunction.h $(UCBDIR)UCBlock.h
	$(CC) -c $*.cpp -o $@ $(SMSINC) $(UCBINC) $(SW)


########################## End of makefile ###################################
