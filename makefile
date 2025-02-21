progs = fft eval_nominal_orbit reward correction shadowing table

OBJS = rk78.o int_rtbp.o prtbp.o fft_module.o utils_module.o cv_module.o \
	   correction_module.o

CFLAGS = -O3 -fPIC #-g -fPIC
LDFLAGS = 
LDLIBS = -lm -lgsl -lgslcblas

all: libhalo.so $(progs)

# Careful!!! Make sure to compile OBJS AND LDLIBS into library, otherwise
# python will not know about GSL functions!

libhalo.so: $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -shared -o $@ $(OBJS) $(LDLIBS)

fft: fft_module.o

eval_nominal_orbit: fft_module.o

reward: rk78.o int_rtbp.o fft_module.o utils_module.o

correction: rk78.o int_rtbp.o prtbp.o utils_module.o cv_module.o correction_module.o

shadowing: rk78.o int_rtbp.o prtbp.o utils_module.o cv_module.o correction_module.o

table: rk78.o int_rtbp.o prtbp.o utils_module.o cv_module.o correction_module.o

clean:
	rm libhalo.so $(OBJS) $(progs)
