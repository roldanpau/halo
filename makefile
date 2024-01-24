halo_functions.so: rk78.c int_rtbp.c
	cc -fPIC -shared -o halo_functions.so rk78.c int_rtbp.c

clean:
	rm halo_functions.so
