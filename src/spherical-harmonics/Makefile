CXXFLAGS = -std=c++0x -g -O3 -fPIC -Wall 
LDFLAGS =

all:	libspherical.so

%.o : %.cc %.h
		$(CXX) $(CXXFLAGS) -c -I. $^

libspherical.so:	SphericalHarmonics.o
		$(CXX) $(CXXFLAGS) -shared -o $@ $(LDFLAGS) $^

lint:	
	python cpplint.py --linelength=120 --counting=detailed ${allfiles}

clean:
	rm *~ libspherical.so *.o

.PHONY:	lint clean
