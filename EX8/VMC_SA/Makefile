CC = g++
CFLAGS = -Wall -O3 --std=c++11

VMC_SA_QM1D.exe : VMC_SA_QM1D.o random.o  
	$(CC) random.o VMC_SA_QM1D.o  -o VMC_SA_QM1D.exe
VMC_SA_QM1D.o : VMC_SA_QM1D.cpp 
	$(CC) -c VMC_SA_QM1D.cpp -o VMC_SA_QM1D.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o VMC_SA_QM1D.exe seed.out

clean_result:
	rm output_Have.dat
	rm output_PsiT2.dat





