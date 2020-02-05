#author: Giorgos Tsalidis 5/2/2020

# define the C/C++ compiler to use
CC = g++
NVCC = nvcc
EXECS = v0 v1 v2 v3 
.PHONY: $(EXECS)
# define command to remove files
RM = rm -f 

all: $(EXECS)

clean:
	$(RM) *.o lib/*.a $(EXECS)

v0:
	$(CC) src/ising-sequential.c validation/main_validate.c -o validate_v0.o  && ./validate_v0.o
v1:
	$(NVCC) src/ising-v1.cu validation/main_validate_cuda.cu -o validate_v1.o && ./validate_v1.o
v2:
	$(NVCC) src/ising-v2.cu validation/main_validate_cuda.cu -o validate_v2.o && ./validate_v2.o
v3:
	$(NVCC) src/ising-v3.cu validation/main_validate_cuda.cu -o validate_v3.o && ./validate_v3.o
