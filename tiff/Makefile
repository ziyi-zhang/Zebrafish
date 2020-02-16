CC = g++
MAINSRCS = tiffReader.cpp 
OTHSRCS =
SRCS = $(MAINSRCS) $(OTHSRCS)
OBJS = $(SRCS:.cpp=.o)
TARGETS = $(MAINSRCS:.cpp=)

CFLAGS = -std=c++11 -Wall
INCPATH = -I./tiff-4.1.0/libtiff
LIBPATH = -L./tiff-4.1.0/libtiff/.libs -ltiff

all: $(TARGETS)

$(TARGETS): $(OBJS)
	$(CC) $(CFLAGS) $(INCPATH) $(OBJS) $(LIBPATH) -o $(TARGETS)

$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -c -o $(OBJS)

memcheck: $(TARGETS)
	valgrind --leak-check=full ./$(TARGETS) $(MEMCHECKARGS) 

clean:
	rm -f $(OBJS) $(TARGETS)

.PHONY: all memcheck clean
