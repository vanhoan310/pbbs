# /home/vanhoan310/server /data/hoan /home/vanhoan310/OptTools
my_dir=/data/hoan/tools
CXX = g++
CFLAGS = -O3 -std=c++11  -fexceptions -DNDEBUG -DIL_STD 
#-pedantic -m64 -O -fPIC -w -fno-strict-aliasing
INCL = -I.
IPATH = -I $(my_dir)/CPLEX_Studio128/cplex/include -I $(my_dir)/CPLEX_Studio128/concert/include -I $(my_dir)/lemon/include -I /home/vanhoan310/OptTools/boost_1_66_0
LPATH = -L $(my_dir)/lemon/lib -lemon -L $(my_dir)/CPLEX_Studio128/cplex/lib/x86-64_linux/static_pic -L $(my_dir)/CPLEX_Studio128/concert/lib/x86-64_linux/static_pic
BINARIES = sc

build: $(BINARIES) 


clean:
	rm -f $(BINARIES) 
	rm -f *.o
	rm -f *~
	rm -f core

# $(BINARIES): test_set_cover_exact_solver.cpp set_cover_exact_solver.cpp set_cover_exact_solver.h 
$(BINARIES): main.cpp util_functions.cpp util_functions.h set_cover_base.cpp set_cover_base.h set_cover_exact_solver.cpp set_cover_exact_solver.h
	$(CXX) $^ -o $@ $(CFLAGS) $(IPATH) $(LPATH) -lpthread -lconcert -lilocplex -lcplex -lm -ldl 

