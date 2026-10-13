module purge
module load gcc/14.2.0 openmpi/5.0.7
module list
cd libpicl
rm ppiclF/source/ppiclf.f
rm ppiclF/source/PPICLF_USER.h
rm ppiclF/source/PPICLF_USER_COMMON.h
cd ppiclF/source
ln -s ../../user_files/PPICLF_USER.h .
ln -s ../../user_files/PPICLF_USER_COMMON.h .
cd ../../
make clean
make
cd .. 
make clean 

rm build_lib/*.f90
rm build_lib/*.d
rm build_lib/*.o

rm build_util/*/*.f90
rm build_util/*/*.d
rm build_util/*/*.o

# Add SPEC=1 only for species transport or the JWL gas model (GASMODEL 7)
make RFLU=1 PICL=1 FOLDER=1 -j16
ls --color=auto
