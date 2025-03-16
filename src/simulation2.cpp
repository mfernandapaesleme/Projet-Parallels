#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>
#include <mpi.h>
#include <fstream>
#include "model.hpp"     // Contient la classe Model et ses méthodes (vegetal_map(), fire_map(), update(), etc.)
#include "display.hpp"   // Pour l'affichage via SDL

using namespace std::string_literals;
using namespace std::chrono_literals;

// Structure pour stocker les paramètres de simulation
struct ParamsType {
    double length {1.0};             // Longueur du terrain (en km)
    unsigned discretization {300u};   // Nombre de cases par direction
    std::array<double,2> wind {0.0, 0.0};  // Vecteur vitesse du vent
    Model::LexicoIndices start {10u,10u};   // Position initiale du foyer (colonne, ligne)
};

// --- Fonctions d'analyse des arguments (identiques à l'exemple précédent) ---
void analyze_arg(int nargs, char* args[], ParamsType& params) {
    if (nargs == 0) return;
    std::string key(args[0]);
    if (key == "-l"s) {
        if (nargs < 2) {
            std::cerr << "Manque une valeur pour la longueur du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.length = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    auto pos = key.find("--longueur=");
    if (pos < key.size()) {
        auto subkey = std::string(key, pos+11);
        params.length = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }
    if (key == "-n"s) {
        if (nargs < 2) {
            std::cerr << "Manque une valeur pour le nombre de cases par direction !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.discretization = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--number_of_cases=");
    if (pos < key.size()) {
        auto subkey = std::string(key, pos+18);
        params.discretization = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }
    if (key == "-w"s) {
        if (nargs < 2) {
            std::cerr << "Manque une paire de valeurs pour le vent !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values(args[1]);
        auto commaPos = values.find(",");
        if (commaPos == std::string::npos) {
            std::cerr << "Deux valeurs séparées par une virgule sont requises pour définir le vent" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.wind[0] = std::stod(values.substr(0, commaPos));
        params.wind[1] = std::stod(values.substr(commaPos+1));
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--wind=");
    if (pos < key.size()) {
        auto subkey = std::string(key, pos+7);
        auto commaPos = subkey.find(",");
        if (commaPos == std::string::npos) {
            std::cerr << "Deux valeurs séparées par une virgule sont requises pour définir le vent" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.wind[0] = std::stod(subkey.substr(0, commaPos));
        params.wind[1] = std::stod(subkey.substr(commaPos+1));
        analyze_arg(nargs-1, &args[1], params);
        return;
    }
    if (key == "-s"s) {
        if (nargs < 2) {
            std::cerr << "Manque une paire de valeurs pour la position initiale du foyer !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values(args[1]);
        auto commaPos = values.find(",");
        if (commaPos == std::string::npos) {
            std::cerr << "Deux valeurs séparées par une virgule sont requises pour définir la position initiale du foyer" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.start.column = std::stoul(values.substr(0, commaPos));
        params.start.row = std::stoul(values.substr(commaPos+1));
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--start=");
    if (pos < key.size()) {
        auto subkey = std::string(key, pos+8);
        auto commaPos = subkey.find(",");
        if (commaPos == std::string::npos) {
            std::cerr << "Deux valeurs séparées par une virgule sont requises pour définir la position initiale du foyer" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.start.column = std::stoul(subkey.substr(0, commaPos));
        params.start.row = std::stoul(subkey.substr(commaPos+1));
        analyze_arg(nargs-1, &args[1], params);
        return;
    }
}

ParamsType parse_arguments(int nargs, char* args[]) {
    if (nargs == 0) return {};
    if (std::string(args[0]) == "--help"s || std::string(args[0]) == "-h"s) {
        std::cout <<
R"RAW(Usage : simulation [option(s)]
  Lance la simulation d'incendie en prenant en compte les options.
  Options :
    -l, --longueur=LONGUEUR     Taille du terrain (en km)
    -n, --number_of_cases=N     Nombre de cases par direction
    -w, --wind=VX,VY            Vecteur vitesse du vent
    -s, --start=COL,ROW         Position initiale du foyer (colonne, ligne)
)RAW";
        exit(EXIT_SUCCESS);
    }
    ParamsType params;
    analyze_arg(nargs, args, params);
    return params;
}

bool check_params(const ParamsType& params) {
    bool flag = true;
    if (params.length <= 0) {
        std::cerr << "[ERREUR FATALE] La longueur du terrain doit être positive et non nulle !" << std::endl;
        flag = false;
    }
    if (params.discretization <= 0) {
        std::cerr << "[ERREUR FATALE] Le nombre de cases par direction doit être positif et non nul !" << std::endl;
        flag = false;
    }
    if (params.start.row >= params.discretization || params.start.column >= params.discretization) {
        std::cerr << "[ERREUR FATALE] Indices incorrects pour la position initiale du foyer" << std::endl;
        flag = false;
    }
    return flag;
}

void display_params(const ParamsType& params) {
    std::cout << "Paramètres de simulation :" << std::endl
              << "\tTaille du terrain : " << params.length << " km" << std::endl
              << "\tNombre de cases par direction : " << params.discretization << std::endl
              << "\tVent : [" << params.wind[0] << ", " << params.wind[1] << "]" << std::endl
              << "\tPosition initiale du foyer (col, ligne) : " 
              << params.start.column << ", " << params.start.row << std::endl;
}

int main(int argc, char* argv[]) {
    // Traitement des arguments de la ligne de commande
    auto params = parse_arguments(argc-1, &argv[1]);
    display_params(params);
    if (!check_params(params)) {
        return EXIT_FAILURE;
    }
    
    // Initialisation de MPI et duplication du communicateur global
    MPI_Init(&argc, &argv);
    MPI_Comm globCom;
    MPI_Comm_dup(MPI_COMM_WORLD, &globCom);
    int nbp;
    MPI_Comm_size(globCom, &nbp);
    int rank;
    MPI_Comm_rank(globCom, &rank);

    int global_size = params.discretization * params.discretization;
    
    // Dans cet exemple, nous supposons que le processus 0 est dédié à l'affichage
    // et le processus 1 (ou tous les autres) font la simulation.
    if (rank == 0) {
        // Processus d'affichage (maître global)
        auto displayer = Displayer::init_instance(params.discretization, params.discretization);
        // Vecteurs pour stocker les cartes globales (feu et végétation)
        std::vector<unsigned char> global_fire_map(global_size);
        std::vector<unsigned char> global_veg_map(global_size);
        
        // Ouverture du fichier pour enregistrer les temps
        std::ofstream file("temps.txt");
        file << "TimeStep Temps_Update Temps_Display Temps_Total\n";
        
        bool loop = true;
        int time_step = 0;
        double start_total = MPI_Wtime();
        
        while(loop) {
            int request = 1;
            // Demande de mise à jour des données au processus de calcul (ici, rang 1)
            MPI_Send(&request, 1, MPI_INT, 1, 0, globCom);
            
            // Réception des cartes globales depuis le processus de calcul
            MPI_Recv(global_fire_map.data(), global_size, MPI_UNSIGNED_CHAR, 1, 0, globCom, MPI_STATUS_IGNORE);
            MPI_Recv(global_veg_map.data(), global_size, MPI_UNSIGNED_CHAR, 1, 1, globCom, MPI_STATUS_IGNORE);
            
            // Réception du temps d'itération (update) depuis le processus de calcul
            double iter_time = 0.0;
            MPI_Recv(&iter_time, 1, MPI_DOUBLE, 1, 2, globCom, MPI_STATUS_IGNORE);
            
            // Mesure du temps d'affichage avec MPI_Wtime()
            double start_display = MPI_Wtime();
            displayer->update(global_veg_map, global_fire_map);
            double end_display = MPI_Wtime();
            
            double display_time = end_display - start_display;
            double total_time = MPI_Wtime() - start_total;
            
            // Enregistrement dans le fichier
            file << time_step << " " 
                 << iter_time << " " 
                 << display_time << " " 
                 << total_time << "\n";
            
            std::cout << "Time step " << time_step << "\n===============" << std::endl;
            
            // Gestion des événements SDL (fermeture de la fenêtre)
            SDL_Event event;
            if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
                int quit_signal = -1;
                MPI_Send(&quit_signal, 1, MPI_INT, 1, 0, globCom);
                loop = false;
            }
            std::this_thread::sleep_for(100ms);
            time_step++;
        }
        file.close();
        SDL_Quit();
    }
    else {
        // Processus de calcul (simulation)
        Model simu(params.length, params.discretization, params.wind, params.start);
        bool continue_simulation = true;
        double total_iter_time = 0.0;
        int iter_count = 0;
        
        while (continue_simulation && simu.update()) {
            double iter_start = MPI_Wtime();
            
            // Vérification non bloquante d'un signal d'arrêt venant du processus d'affichage
            int flag = 0;
            MPI_Status status;
            MPI_Iprobe(0, 0, globCom, &flag, &status);
            if (flag) {
                int signal;
                MPI_Recv(&signal, 1, MPI_INT, 0, 0, globCom, MPI_STATUS_IGNORE);
                if (signal == -1) {
                    continue_simulation = false;
                    break;
                }
            }
            
            // Envoi des cartes de simulation (feu et végétation) au processus d'affichage
            MPI_Send(simu.fire_map().data(), global_size, MPI_UNSIGNED_CHAR, 0, 0, globCom);
            MPI_Send(simu.vegetal_map().data(), global_size, MPI_UNSIGNED_CHAR, 0, 1, globCom);
            
            double iter_end = MPI_Wtime();
            double iter_time = iter_end - iter_start;
            total_iter_time += iter_time;
            iter_count++;
            
            // Envoi du temps d'itération au processus d'affichage (tag 2)
            MPI_Send(&iter_time, 1, MPI_DOUBLE, 0, 2, globCom);
            
            std::this_thread::sleep_for(100ms);
        }
        
        double average_iter_time = (iter_count > 0) ? (total_iter_time / iter_count) : 0.0;
        std::cout << "Temps moyen par itération (calcul + communication) : " 
                  << average_iter_time << " secondes" << std::endl;
    }
    
    MPI_Finalize();
    return EXIT_SUCCESS;
}
