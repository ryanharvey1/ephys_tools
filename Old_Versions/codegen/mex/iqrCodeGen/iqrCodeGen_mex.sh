MATLAB="/Applications/MATLAB_R2015a.app"
Arch=maci64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/Users/ryanharvey/.matlab/R2015a"
OPTSFILE_NAME="./setEnv.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for iqrCodeGen" > iqrCodeGen_mex.mki
echo "CC=$CC" >> iqrCodeGen_mex.mki
echo "CFLAGS=$CFLAGS" >> iqrCodeGen_mex.mki
echo "CLIBS=$CLIBS" >> iqrCodeGen_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> iqrCodeGen_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> iqrCodeGen_mex.mki
echo "CXX=$CXX" >> iqrCodeGen_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> iqrCodeGen_mex.mki
echo "CXXLIBS=$CXXLIBS" >> iqrCodeGen_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> iqrCodeGen_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> iqrCodeGen_mex.mki
echo "LD=$LD" >> iqrCodeGen_mex.mki
echo "LDFLAGS=$LDFLAGS" >> iqrCodeGen_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> iqrCodeGen_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> iqrCodeGen_mex.mki
echo "Arch=$Arch" >> iqrCodeGen_mex.mki
echo OMPFLAGS= >> iqrCodeGen_mex.mki
echo OMPLINKFLAGS= >> iqrCodeGen_mex.mki
echo "EMC_COMPILER=Xcode with Clang" >> iqrCodeGen_mex.mki
echo "EMC_CONFIG=optim" >> iqrCodeGen_mex.mki
"/Applications/MATLAB_R2015a.app/bin/maci64/gmake" -B -f iqrCodeGen_mex.mk
