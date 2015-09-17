CC = g++ -std=gnu++11 -Wall
LIB = -lntl
HDRS = -I gnfs
OBJS = gnfs/NTL_extension.o gnfs/block_lanczos.cpp gnfs/polynomial_selection.o gnfs/sieve.o gnfs/linear_algebra.o gnfs/square_root.o

default: $(OBJS)
	$(CC) $(HDRS) gnfs/prime_gen.cpp -o primegen $(LIB)
	$(CC) $(HDRS) gnfs/n_gen.cpp -o ngen $(LIB)
	$(CC) $(HDRS) $(OBJS) gnfs/gnfs.cpp -lntl -o factor_gnfs $(LIB)


%.o: %.cpp %.hpp
	$(CC) $(HDRS) -c $< -o $@

clean:
	rm -f gnfs/*.o
	rm -f primegen
	rm -f factor_gnfs
	rm -f ngen
