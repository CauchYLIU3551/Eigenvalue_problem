# by Chengyu Liu 16th June

include ~/AFEPack/Make.global_options

source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))
LDFLAGS += -L~/AFEPack/library/lib -lAFEPack

all : main

%.o : %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main : $(object)
	$(CXX) -o $@ $(object) $(LDFLAGS) $(LIBS)

clean :
	-rm -rf $(object)
	-rm -rf main
	-rm -f *.[nes]
	-rm -f *.dx

.PHONY : default clean
