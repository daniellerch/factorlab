
CC = g++ -std=gnu++11 -Wall
LIB = -lntl
HDRS = -I src

OBJS =   NTL_extension.o \
         \
         polynomial_evaluation.o \
         pollard_rho.o \
         pollard_rho_ppa.o \
         pollard_pm1.o \
         quadratic_sieve.o \
         \
         block_lanczos.o \
         polynomial_selection.o \
         sieve.o \
         linear_algebra.o \
         square_root.o \
    
BUILD = $(addprefix build/,$(OBJS))    

default: $(addprefix src/,$(OBJS))
	$(CC) $(HDRS) $(BUILD) src/quadratic_sieve_main.cpp -o bin/factor_quadratic_sieve $(LIB)
	$(CC) $(HDRS) $(BUILD) src/polynomial_evaluation_main.cpp -o bin/factor_polynomial_evaluation $(LIB)
	$(CC) $(HDRS) $(BUILD) src/pollard_pm1_main.cpp -o bin/factor_pollard_pm1 $(LIB)
	$(CC) $(HDRS) $(BUILD) src/pollard_rho_main.cpp -o bin/factor_pollard_rho $(LIB)
	$(CC) $(HDRS) $(BUILD) src/n_gen_main.cpp -o bin/n_gen $(LIB)
	$(CC) $(HDRS) $(BUILD) src/prime_gen_main.cpp -o bin/prime_gen $(LIB)
	#$(CC) $(HDRS) $(BUILD)) src/pollard_rho_ppa_main.cpp -o bin/factor_pollard_rho_ppa $(LIB)
	$(CC) $(HDRS) $(BUILD) src/gnfs.cpp -lntl -o bin/factor_gnfs $(LIB)


%.o: %.cpp
	$(CC) $(HDRS) -c $(subst build,src,$<) -o $(subst src,build,$@)

clean:
	rm -f build/*.o bin/*

