CPPFLAGS = 

walk: dynein_walk.cpp
	g++ dynein_walk.cpp -o dynein_walk $(CPPFLAGS)

clean:
	rm -f *.o
	rm -f dynein_walk

#	rm -f prevents rm from complaining when there is no such file to delete
