setsm : Function.o 
	gcc -O3 -ffast-math -funroll-loops -fopenmp -I/usr/lib/x86_64-redhat-linux5E/include -o setsm Function.o -lm -ltiff

Function.o : Typedefine.h Function.h Function.c
	gcc -O3 -ffast-math -funroll-loops -fopenmp -I/usr/lib/x86_64-redhat-linux5E/include -c -static Typedefine.h Function.h Function.c -lm -ltiff 