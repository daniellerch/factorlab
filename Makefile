
CC = g++ -std=gnu++11 -Wall
LIB = -lntl
HDRS = -I src

OBJS =   src/NTL_extension.o \
         src/pollard_rho.o \
         src/pollard_pm1.o \
         src/block_lanczos.cpp \
         src/polynomial_selection.o \
         src/sieve.o \
         src/linear_algebra.o \
         src/square_root.o

default: $(OBJS)
	$(CC) $(HDRS) $(OBJS) src/pollard_pm1_main.cpp -o factor_pollard_pm1 $(LIB)
	$(CC) $(HDRS) $(OBJS) src/pollard_rho_main.cpp -o factor_pollard_rho $(LIB)
	$(CC) $(HDRS) src/prime_gen.cpp -o primegen $(LIB)
	$(CC) $(HDRS) src/n_gen.cpp -o ngen $(LIB)
	$(CC) $(HDRS) $(OBJS) src/gnfs.cpp -lntl -o factor_gnfs $(LIB)


%.o: %.cpp
	$(CC) $(HDRS) -c $< -o $@

clean:
	rm -f src/*.o
	rm -f primegen
	rm -f factor_gnfs
	rm -f ngen
