progs = fft reward correction

OBJS = rk78.o int_rtbp.o fft_module.o utils_module.o cv_module.o

CFLAGS = -g -fPIC #-O3
LDFLAGS = 
LDLIBS = -lm -lgsl -lgslcblas

all: libhalo.so $(progs)

# Careful!!! Make sure to compile OBJS AND LDLIBS into library, otherwise
# python will not know about GSL functions!

libhalo.so: $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -shared -o $@ $(OBJS) $(LDLIBS)

fft: fft_module.o

reward: rk78.o int_rtbp.o fft_module.o utils_module.o

correction: rk78.o int_rtbp.o utils_module.o cv_module.o

clean:
	rm libhalo.so $(OBJS) $(progs)
