all:
	g++ -Wall -g -o run.exe ./main.cxx ./EffStudy.cxx `root-config --libs --glibs --cflags` -lTreePlayer

clean:
	rm run.exe -f *.o *.so *.d

