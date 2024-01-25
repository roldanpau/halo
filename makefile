progs = fft

OBJS = rk78.o int_rtbp.o fft_module.o

CFLAGS = -g -fPIC #-O3
LDFLAGS = 
LDLIBS = -lm -lgsl -lgslcblas

all: libhalo.so $(progs)

# Careful!!! Make sure to compile OBJS AND LDLIBS into library, otherwise
# python will not know about GSL functions!

libhalo.so: $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -shared -o $@ $(OBJS) $(LDLIBS)

fft: fft_module.o

clean:
	rm libhalo.so $(OBJS) $(progs)
