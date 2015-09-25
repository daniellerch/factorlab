
CC = g++ -std=gnu++11 -Wall
LIB = -lntl
HDRS = -I src -I src/gnfs

OBJS =   src/ntl_ext/NTL_extension.o \
         \
         src/poly_eval/polynomial_evaluation.o \
         src/rho/pollard_rho.o \
         src/rho/ec_rho.o \
         src/pm1/pollard_pm1.o \
         src/ecm/ecm.o \
         \
         src/curves/EC_p.o \
         src/curves/pair_ZZ_long.o \
         src/matrix/SGauss.o \
         src/matrix/mat_long.o \
         src/matrix/svector.o \
         src/matrix/smatrix.o \
         src/matrix/block_lanczos.o \
         \
         src/qs/quadratic_sieve.o \
         \
         src/gnfs/polynomial_selection.o \
         src/gnfs/sieve.o \
         src/gnfs/linear_algebra.o \
         src/gnfs/square_root.o \
    

default: $(OBJS)
	# -------------------------
	# Factorization Algorithms 
	# -------------------------
	$(CC) $(HDRS) $(OBJS) src/ecm/ecm_main.cpp -o bin/factor_ecm $(LIB)
	$(CC) $(HDRS) $(OBJS) src/qs/quadratic_sieve_main.cpp -o bin/factor_quadratic_sieve $(LIB)
	$(CC) $(HDRS) $(OBJS) src/poly_eval/polynomial_evaluation_main.cpp -o bin/factor_polynomial_evaluation $(LIB)
	$(CC) $(HDRS) $(OBJS) src/pm1/pollard_pm1_main.cpp -o bin/factor_pollard_pm1 $(LIB)
	$(CC) $(HDRS) $(OBJS) src/rho/pollard_rho_main.cpp -o bin/factor_pollard_rho $(LIB)
	$(CC) $(HDRS) $(OBJS) src/gnfs/gnfs_main.cpp -o bin/factor_gnfs $(LIB)
	# -------------------------
	# Tools
	# -------------------------
	$(CC) $(HDRS) $(OBJS) src/tools/ec_rho_main.cpp -o bin/factor_ec_rho $(LIB)
	$(CC) $(HDRS) $(OBJS) src/tools/n_gen_main.cpp -o bin/n_gen $(LIB)
	$(CC) $(HDRS) $(OBJS) src/tools/prime_gen_main.cpp -o bin/prime_gen $(LIB)

%.o: %.cpp
	$(CC) $(HDRS) -c $< -o $@

clean:
	rm -f src/*.o src/*/*.o bin/*

