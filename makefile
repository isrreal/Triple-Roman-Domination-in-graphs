CPPFLAGS=-std=c++17 -Wall -Wextra -O3
IPATH=-Iinc/
SRC_FOLDER=src/
OBJ_FOLDER=obj/

SOURCES= $(SRC_FOLDER)main.cpp $(SRC_FOLDER)GeneticAlgorithm.cpp $(SRC_FOLDER)Chromosome.cpp \
         $(SRC_FOLDER)Graph.cpp $(SRC_FOLDER)TripleRomanDomination.cpp $(SRC_FOLDER)AntColonyOptimization.cpp \
         $(SRC_FOLDER)util_functions.cpp

OBJECTS= $(SOURCES:$(SRC_FOLDER)%.cpp=$(OBJ_FOLDER)%.gch)

all: create_obj_dir app

app: $(OBJECTS)
	g++ $(OBJECTS) -o app

$(OBJ_FOLDER)%.gch: $(SRC_FOLDER)%.cpp
	g++ $(CPPFLAGS) $(IPATH) -c $< -o $@

create_obj_dir:
	mkdir -p $(OBJ_FOLDER)
	
.PHONY: clean_cache

clean_cache:
	ccache --clear

clean:
	rm -rf $(OBJ_FOLDER) app

