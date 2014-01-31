CPPFLAGS = 

walk: dynein_walk.cpp
	g++ dynein_walk.cpp -o walk $(CPPFLAGS)

clean:
	rm -f *.o
	rm -f walk
	rm data.txt

#	rm -f prevents rm from complaining when there is no such file to delete
