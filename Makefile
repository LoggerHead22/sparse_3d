CXX := g++
CXXFLAGS := -std=c++17 -ffast-math -march=native  -mfpmath=sse -O3 #-W -Wuninitialized -Wall -Wunused -fsanitize=address -mfpmath=sse -fstack-protector-all -g   -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
RM := del

a.out: main.o mls_apr.o 
	$(CXX) $^ -o $@
	
%.o: %.cpp *.hpp
	$(CXX) -c $(CXXFLAGS) $< -o $@
	
clean:
	$(RM) *.o *.exe