CC:=icc
FC:=ifort
#Generic AMD
ifeq ($(CPU),amd)
    OPTFLAGS:=-ipo -O3 -msse3 -no-prec-div
#Generic Intel
else ifeq ($(CPU),intel)
    OPTFLAGS:=-ipo -O3 -no-prec-div -xSSE4.2
# Host
else
    OPTFLAGS:=-ipo -O3 -xHost
endif
OPTFLAGS += -heap-arrays
DBGFLAGS:=-O0 -g 
FDBG:=-fpe0 -traceback -check all -ftrapuv -warn unused
FFLAGS:=
#CPPFLAGS += -I/opt/local/include/ufsparse
CPPFLAGS += -I/usr/include/suitesparse

LAPACK := -mkl=parallel
#LIBS += -lcholmod -lamd -lcamd -lcolamd -lccolamd
