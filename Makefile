GENERIC_FLAGS = -Wall -O3 
UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
LIBS =  -lm -lpthread -lfftw3 -lfftw3_threads -framework OpenGL -framework GLUT -lobjc
FLAGS = -DMAC $(GENERIC_FLAGS)
else
LIBS =  -lm -lpthread -lfftw3 -lfftw3_threads -lX11 -lXext -lGL -lGLU -lglut
FLAGS = $(GENERIC_FLAGS)
endif


wavy: particle.o main.o 
	gcc -o wavy particle.o main.o $(LIBS) -L/opt/homebrew/Cellar/fftw/3.3.11/lib

%.o: %.c
	gcc $(FLAGS) -o $(@) -c $*.c -I/opt/homebrew/Cellar/fftw/3.3.11/include/ 

clean :
	rm -f *.o core wavy 

