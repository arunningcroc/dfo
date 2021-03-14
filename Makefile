OBJS = types.o basis.o hamilton.o solver.o energies.o main.o
OBJS2 = types.o basis.o energies.o optimize.o minim.o
OBJS3 = types.o basis.o hamilton.o energies.o solver.o optimize.o test.o
FC = gfortran

#FFLAGS = -Og -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation \
#         -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter \
#         -fwhole-file -fcheck=all -std=f2008 -pedantic -fbacktrace \
#         -finit-real=nan -Wextra -Wconversion -pg
#-O3 -funroll-loops
FFLAGS = -march=native -mtune=native -O2 -fopenmp
#FFLAGS = -O0

LIBS = -lblas -llapack

PROGRAM = main
MINI = minim
TESTI = test

all: $(PROGRAM)
minimizer: $(MINI)
testaa: $(TESTI)

$(TESTI): $(OBJS3)
	$(FC) $(FFLAGS) $(OBJS3) -o $(TESTI) $(LIBS)
$(MINI): $(OBJS2)
	$(FC) $(FFLAGS) $(OBJS2) -o $(MINI) $(LIBS)
$(PROGRAM): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(PROGRAM) $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< $(LIBS)
clean:
	$(RM) *.mod *.o *.dt fort.*