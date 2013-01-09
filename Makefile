LIBS =  -lm -lpthread -lfftw3 -lfftw3_threads -lX11 -lXext -lGL -lGLU -lglut
OPT = -O4
DEBUG = -g -DDEBUG
#FLAGS = -Wall $(DEBUG) 
FLAGS = -Wall $(OPT)

wavy: main.o 
	gcc -o wavy main.o $(LIBS)

%.o: %.c
	$(CC) $(FLAGS) -o $(@) -c $*.c 

clean :
	rm -f *.o core wavy 

