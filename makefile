CPPFLAGS=-std=c++17 
IPATH=-Iinc/
SRC=src/
OBJ=obj/

all: create_obj_dir app

app: main.gch GeneticAlgorithm.gch Chromosome.gch Graph.gch TripleRomanDomination.gch AntColonyOptimization.gch
	g++ $(OBJ)main.gch $(OBJ)GeneticAlgorithm.gch $(OBJ)Chromosome.gch $(OBJ)Graph.gch $(OBJ)TripleRomanDomination.gch $(OBJ)AntColonyOptimization.gch -o app

main.gch: $(SRC)main.cpp
	g++ $(CPPFLAGS) $(IPATH) -c $(SRC)main.cpp -o $(OBJ)main.gch

GeneticAlgorithm.gch: $(SRC)GeneticAlgorithm.cpp                         
	g++ $(CPPFLAGS) $(IPATH) -c $(SRC)GeneticAlgorithm.cpp -o $(OBJ)GeneticAlgorithm.gch

AntColonyOptimization.gch: $(SRC)AntColonyOptimization.cpp                         
	g++ $(CPPFLAGS) $(IPATH) -c $(SRC)AntColonyOptimization.cpp -o $(OBJ)AntColonyOptimization.gch

TripleRomanDomination.gch: $(SRC)TripleRomanDomination.cpp 
	g++ $(CPPFLAGS) $(IPATH) -c $(SRC)TripleRomanDomination.cpp -o $(OBJ)TripleRomanDomination.gch

Chromosome.gch: $(SRC)Chromosome.cpp 
	g++ $(CPPFLAGS) $(IPATH) -c $(SRC)Chromosome.cpp -o $(OBJ)Chromosome.gch
	
Graph.gch: $(SRC)Graph.cpp 
	g++ $(CPPFLAGS) $(IPATH) -c $(SRC)Graph.cpp -o $(OBJ)Graph.gch

create_obj_dir:
	mkdir -p $(OBJ)

clean:
	rm -rf $(OBJ)*.gch app

