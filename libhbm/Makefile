all: libhbm.so

libhbm.so: libhbm.cpp
	icpc -shared libhbm.cpp -fPIC -ltbbmalloc -lnuma -o libhbm.so -ldl

clean:
	rm libhbm.so
