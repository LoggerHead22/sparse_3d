CXX := g++
CXXFLAGS := -std=c++17 -Wall   #-ffast-math -march=native  -mfpmath=sse 
RM := del

a.out: main.o mls_apr.o 
	$(CXX) $^ -o $@
	
%.o: %.cpp *.hpp
	$(CXX) -c $(CXXFLAGS) $< -o $@
	
clean:
	$(RM) *.o *.exe