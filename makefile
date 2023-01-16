all:
	g++ main.cpp Izing.cpp -pthread -o model.exe
	./model.exe