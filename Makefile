
CC = g++ -std=gnu++11 -Wall
LIB = -lntl
HDRS = -I src

OBJS =   src/NTL_extension.o \
         \
         src/polynomial_evaluation.o \
         src/pollard_rho.o \
         src/pollard_pm1.o \
         src/quadratic_sieve.o \
         src/ecm.o \
         \
         src/ec/EC_p.o \
         src/ec/pair_ZZ_long.o \
         src/matrix/SGauss.o \
         src/matrix/mat_long.o \
         src/matrix/svector.o \
         src/matrix/smatrix.o \
         src/matrix/block_lanczos.o \
         src/polynomial_selection.o \
         src/sieve.o \
         src/linear_algebra.o \
         src/square_root.o \
    

default: $(OBJS)
	# -------------------------
	# Factorization Algorithms 
	# -------------------------
	$(CC) $(HDRS) $(OBJS) src/ecm_main.cpp -o bin/factor_ecm $(LIB)
	$(CC) $(HDRS) $(OBJS) src/quadratic_sieve_main.cpp -o bin/factor_quadratic_sieve $(LIB)
	$(CC) $(HDRS) $(OBJS) src/polynomial_evaluation_main.cpp -o bin/factor_polynomial_evaluation $(LIB)
	$(CC) $(HDRS) $(OBJS) src/pollard_pm1_main.cpp -o bin/factor_pollard_pm1 $(LIB)
	$(CC) $(HDRS) $(OBJS) src/pollard_rho_main.cpp -o bin/factor_pollard_rho $(LIB)
	$(CC) $(HDRS) $(OBJS) src/gnfs_main.cpp -o bin/factor_gnfs $(LIB)
	# -------------------------
	# Tools
	# -------------------------
	$(CC) $(HDRS) $(OBJS) src/n_gen_main.cpp -o bin/n_gen $(LIB)
	$(CC) $(HDRS) $(OBJS) src/prime_gen_main.cpp -o bin/prime_gen $(LIB)

%.o: %.cpp
	$(CC) $(HDRS) -c $< -o $@

clean:
	rm -f src/*.o src/*/*.o bin/*

