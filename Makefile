
CC = g++ -std=gnu++11 -Wall
LIB = -lntl
HDRS = -I src -I src/gnfs

OBJS_GNFS =   src/ntl_ext/NTL_extension.o \
         \
         src/matrix/SGauss.o \
         src/matrix/mat_long.o \
         src/matrix/svector.o \
         src/matrix/smatrix.o \
         src/matrix/block_lanczos.o \
         \
         src/gnfs/polynomial_selection.o \
         src/gnfs/sieve.o \
         src/gnfs/linear_algebra.o \
         src/gnfs/square_root.o \
          
   
OBJS_POE =  src/poly_eval/polynomial_evaluation.o 

OBJS_RHO =  src/rho/pollard_rho.o 

OBJS_PM1 =  src/pm1/pollard_pm1.o 

OBJS_ECM =  src/ecm/ecm.o \
            src/curves/EC_p.o \
            src/curves/pair_ZZ_long.o \


OBJS_ECR =  src/ecm/ecm.o \
            src/rho/ec_rho.o \
            src/curves/EC_p.o \
            src/curves/pair_ZZ_long.o \

OBJS_QS = src/qs/quadratic_sieve.o             

OBJS = $(OBJS_GNFS) $(OBJS_QS) $(OBJS_POE) $(OBJS_RHO) $(OBJS_PM1) $(OBJS_ECM) $(OBJS_ECR)

default: $(OBJS)
	# -------------------------
	# Factorization Algorithms 
	# -------------------------
	$(CC) $(HDRS) $(OBJS_ECM) src/ecm/ecm_main.cpp -o bin/factor_ecm $(LIB)
	$(CC) $(HDRS) $(OBJS_QS) src/qs/quadratic_sieve_main.cpp -o bin/factor_quadratic_sieve $(LIB)
	$(CC) $(HDRS) $(OBJS_POE) src/poly_eval/polynomial_evaluation_main.cpp -o bin/factor_polynomial_evaluation $(LIB)
	$(CC) $(HDRS) $(OBJS_PM1) src/pm1/pollard_pm1_main.cpp -o bin/factor_pollard_pm1 $(LIB)
	$(CC) $(HDRS) $(OBJS_RHO) src/rho/pollard_rho_main.cpp -o bin/factor_pollard_rho $(LIB)
	$(CC) $(HDRS) $(OBJS_GNFS) src/gnfs/gnfs_main.cpp -o bin/factor_gnfs $(LIB)
	# -------------------------
	# Tools
	# -------------------------
	$(CC) $(HDRS) $(OBJS_ECR) src/tools/ec_rho_main.cpp -o bin/factor_ec_rho $(LIB)
	$(CC) $(HDRS) src/tools/n_gen_main.cpp -o bin/n_gen $(LIB)
	$(CC) $(HDRS) src/tools/prime_gen_main.cpp -o bin/prime_gen $(LIB)

%.o: %.cpp
	$(CC) $(HDRS) -c $< -o $@

clean:
	rm -f src/*.o src/*/*.o bin/*

