include Make_linux.inc
#include Make_msys2.inc
#include Make_osx.inc

# Compilateurs
CXX = g++
MPICXX = mpicxx

# Flags de compilation
CXXFLAGS = -std=c++17 -Wall -O3
CXXFLAGS_SEQ = $(CXXFLAGS) -fopenmp
CXXFLAGS_PAR = $(CXXFLAGS) -fopenmp
CXXFLAGS_MPI_PAR = $(CXXFLAGS) -fopenmp

# Flags de linkage
LDFLAGS_SEQ = -lSDL2 -fopenmp
LDFLAGS_PAR = -lSDL2 -fopenmp
LDFLAGS_MPI_PAR = -lSDL2 -fopenmp

# Executables
SEQ_EXEC = simulation_seq.exe
PAR_EXEC = simulation_par.exe

# Objets pour exécutable séquentiel
SEQ_OBJECTS = simulation_seq.o model.o display.o

# Objets pour exécutable parallèle
PAR_OBJECTS = simulation_par.o model1.o display.o

# Cible par défaut : compiler les deux exécutables
all: $(SEQ_EXEC) $(PAR_EXEC)

# Règles pour l'exécutable séquentiel
$(SEQ_EXEC): $(SEQ_OBJECTS) 
		$(CXX) $(CXXFLAGS_SEQ) $^ -o $@ $(LDFLAGS_SEQ)

simulation_seq.o: simulation.cpp model.cpp model.hpp display.hpp
		$(CXX) $(CXXFLAGS_SEQ) -c simulation.cpp -o simulation_seq.o

model.o: model.cpp model.hpp
		$(CXX) $(CXXFLAGS_SEQ) -c model.cpp -o model.o

# Règles pour l'exécutable parallèle
$(PAR_EXEC): $(PAR_OBJECTS)
		$(MPICXX) $(CXXFLAGS_MPI_PAR) $^ -o $@ $(LDFLAGS_MPI_PAR)

simulation_par.o: simulation3.cpp model1.cpp model.hpp display.hpp
		$(MPICXX) $(CXXFLAGS_PAR) -c simulation3.cpp -o simulation_par.o

model1.o: model1.cpp model.hpp
		$(MPICXX) $(CXXFLAGS_PAR) -c model1.cpp -o model1.o

# Objets partagés (display)
display.o: display.cpp display.hpp
		$(CXX) $(CXXFLAGS_SEQ) -c display.cpp -o display.o

# Cible pour nettoyer les fichiers compilés
clean:
		rm -f *.o $(SEQ_EXEC) $(PAR_EXEC) *.log temps.txt

# Cible d'aide
help:
		@echo "Cibles disponibles :"
		@echo "  all                - Compiler les exécutables séquentiel et parallèle"
		@echo "  clean              - Supprimer les fichiers compilés"
		@echo "  run_seq            - Exécuter la simulation séquentielle"
		@echo "  run_parallel_*     - Exécuter les simulations parallèles avec différents nombres de threads"
		@echo "  run_all_parallel_tests - Exécuter tous les tests parallèles"
		@echo "  analyze            - Analyser les résultats des tests"

# Cibles pour exécuter les simulations

# Exécuteur séquentiel
run_seq:
		./$(SEQ_EXEC) -n 300 -l 1.0 -w 0.0,0.0 -s 150,150 -t 1 > sequential.log 2>&1

# Exécuteurs parallèles avec différents nombres de threads
run_parallel_2:
		mpiexec --bind-to none -n 2 ./$(PAR_EXEC) -n 300 -l 1.0 -w 0.0,0.0 -s 150,150 -t 2 > parallel_2.log 2>&1

run_parallel_4:
		mpiexec --bind-to none --oversubscribe -n 4 ./$(PAR_EXEC) -n 300 -l 1.0 -w 0.0,0.0 -s 150,150 -t 4 > parallel_4.log 2>&1

run_parallel_8:
		mpiexec --bind-to none --oversubscribe -n 8 ./$(PAR_EXEC) -n 300 -l 1.0 -w 0.0,0.0 -s 150,150 -t 8 > parallel_8.log 2>&1

run_parallel_16:
		mpiexec --bind-to none --oversubscribe -n 16 ./$(PAR_EXEC) -n 300 -l 1.0 -w 0.0,0.0 -s 150,150 -t 16 > parallel_16.log 2>&1

# Cible pour exécuter tous les tests parallèles
run_all_parallel_tests: run_parallel_2 run_parallel_4 run_parallel_8 run_parallel_16

# Cible pour exécuter le script d'analyse après les tests
analyze:
		./analyze_results.sh

# Cible pour compiler et exécuter tout avec analyse
all_with_analysis: all run_all_parallel_tests analyze