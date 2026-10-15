FC               = mpif77 # Valid MPI Fortran compiler
CC               = mpicc  # Valid MPI C compiler
FFLAGS          += -cpp -w -O3 -std=legacy -g -fmax-errors=10 -fallow-argument-mismatch -fcheck=all -fbacktrace
CPFLAGS         += -E
 	           
INSTALL_LOCATION = .

# Fortran module files (.mod) of ppiclF go here, not in libpicl/
MOD_DIR = $(INSTALL_LOCATION)/build

####################
# DO NOT TOUCH BELOW
####################

# PPICLF LIBRARY
SOURCE_ROOT_PPICLF=$(INSTALL_LOCATION)/ppiclF/source
SOURCE_ROOT_USER=$(INSTALL_LOCATION)/user_files
PPICLF_IFLAGS+=-I$(SOURCE_ROOT_PPICLF)

# GSLIB LIBRARY
SOURCE_ROOT_GSLIB=$(INSTALL_LOCATION)/ppiclF/3rd_party/gslib
SOURCE_ROOT_GSLIB_OBJ=$(SOURCE_ROOT_GSLIB)/gslib/src
GSLIB_IFLAGS+=-I$(SOURCE_ROOT_GSLIB)/include

SRC = $(SOURCE_ROOT_PPICLF)/ppiclf_unsteady.f \
      $(SOURCE_ROOT_PPICLF)/ppiclf_data.f \
      $(SOURCE_ROOT_USER)/ppiclf_user.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_subbinMap.f       \
      $(SOURCE_ROOT_USER)/ppiclf_user_QuasiSteady.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_Fluctuations.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_PseudoTurbulence.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_AddedMass.f     \
      $(SOURCE_ROOT_USER)/Added_Mass/ppiclf_user_AddedMass_functions.f   \
      $(SOURCE_ROOT_USER)/Added_Mass/ppiclf_user_AddedMass_subroutines.f   \
      $(SOURCE_ROOT_USER)/ppiclf_user_HeatTransfer.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_Torque.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_Lift.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_ViscousUnsteady.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_EvalNearestNeighbor.f     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_comm.f     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_op.f       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_io.f       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_solve.f    \
      $(SOURCE_ROOT_PPICLF)/ppiclf_wall.f     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_mxm.f      
OBJ = $(SOURCE_ROOT_PPICLF)/ppiclf_unsteady.o \
      $(SOURCE_ROOT_PPICLF)/ppiclf_data.o \
      $(SOURCE_ROOT_USER)/ppiclf_user.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_subbinMap.o       \
      $(SOURCE_ROOT_USER)/ppiclf_user_QuasiSteady.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_Fluctuations.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_PseudoTurbulence.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_AddedMass.o     \
      $(SOURCE_ROOT_USER)/Added_Mass/ppiclf_user_AddedMass_functions.o   \
      $(SOURCE_ROOT_USER)/Added_Mass/ppiclf_user_AddedMass_subroutines.o   \
      $(SOURCE_ROOT_USER)/ppiclf_user_HeatTransfer.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_Torque.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_Lift.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_ViscousUnsteady.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_EvalNearestNeighbor.o     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_comm.o     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_op.o       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_io.o       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_solve.o    \
      $(SOURCE_ROOT_PPICLF)/ppiclf_wall.o     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_mxm.o
PRE = $(SOURCE_ROOT_PPICLF)/PPICLF_SOLN.h     \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GRID.h     \
      $(SOURCE_ROOT_PPICLF)/PPICLF_OPT.h      \
      $(SOURCE_ROOT_PPICLF)/PPICLF_PARALLEL.h \
      $(SOURCE_ROOT_PPICLF)/PPICLF_USER_COMMON.h \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GEOM.h
POS = $(SOURCE_ROOT_PPICLF)/PPICLF
      

# Make commands
default: makeThird preProcess getObjs libObjs 

libObjs: $(SOURCE_ROOT_PPICLF)/ppiclf.o
	@ar crv $(SOURCE_ROOT_PPICLF)/libppiclF.a $(SOURCE_ROOT_PPICLF)/ppiclf.o $(SOURCE_ROOT_GSLIB_OBJ)/*.o
	@echo "                       "
	@echo "***********************"
	@echo "*** LIBRARY SUCCESS ***"
	@echo "***********************"
	@echo "                       "

getObjs: $(SOURCE_ROOT_PPICLF)/ppiclf.f
	mkdir -p $(MOD_DIR)
	$(FC) $(FFLAGS) -J $(MOD_DIR) -c $(SOURCE_ROOT_PPICLF)/ppiclf.f $(GSLIB_IFLAGS) 
	mv *.o $(SOURCE_ROOT_PPICLF)
	@echo "                              "
	@echo "******************************"
	@echo "*** LIBRARY OBJECT SUCCESS ***"
	@echo "******************************"
	@echo "                              "

preProcess: $(PRE)
	$(CC) $(CPFLAGS) $(PRE) > $(POS)_tmp
	grep "^[^#]" $(POS)_tmp > $(POS)
	rm $(POS)_tmp
	echo '#include "PPICLF_USER.h"' > $(SOURCE_ROOT_PPICLF)/ppiclf.f
	echo '#include "PPICLF_STD.h"' >> $(SOURCE_ROOT_PPICLF)/ppiclf.f
	cat $(SRC)  >> $(SOURCE_ROOT_PPICLF)/ppiclf.f
	@echo "                       "
	@echo "********************"
	@echo "*** PREPROCESSED ***"
	@echo "********************"
	@echo "                       "
	

makeThird: $(SOURCE_ROOT_GSLIB)/install
	cd $(SOURCE_ROOT_GSLIB); \
	./install $(CC) $(FC)
	@echo "                       "
	@echo "***********************"
	@echo "*** INSTALLED GSLIB ***"
	@echo "***********************"
	@echo "                       "

clean:
	 rm -rf $(MOD_DIR)
	 rm -r $(SOURCE_ROOT_PPICLF)/*.o                        
	 rm -r $(SOURCE_ROOT_PPICLF)/ppiclf.f
	 rm -r $(SOURCE_ROOT_PPICLF)/PPICLF
	 rm -r $(SOURCE_ROOT_PPICLF)/libppiclF.a                
	 rm -r $(SOURCE_ROOT_PPICLF)/PPICLF_tmp/*                
	 rm -r $(SOURCE_ROOT_GSLIB)/gslib                       
	 rm -r $(SOURCE_ROOT_GSLIB)/lib                         
	 rm -r $(SOURCE_ROOT_GSLIB)/include                     
	 rm -r $(SOURCE_ROOT_GSLIB)/*.tar.gz                    

