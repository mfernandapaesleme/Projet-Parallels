#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>
#include <mpi.h>
#include "model.hpp"
#include "display.hpp"

using namespace std::string_literals;
using namespace std::chrono_literals;

// Structure pour stocker les paramètres de simulation
struct ParamsType {
    double length {1.0};             // Longueur du terrain (en km)
    unsigned discretization {300u};   // Nombre de cases par direction
    std::array<double,2> wind {0.0, 0.0};  // Vecteur vitesse du vent
    Model::LexicoIndices start {10u,10u};   // Position initiale du foyer (colonne, ligne)
};

// Fonction récursive d'analyse des arguments de la ligne de commande
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

// Fonction de parsing des arguments
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

// Vérification de la validité des paramètres
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

// Affichage des paramètres de simulation
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
    
    // Initialisation de MPI
    MPI_Init(&argc, &argv);
    MPI_Comm globCom;
    MPI_Comm_dup(MPI_COMM_WORLD, &globCom);
    int nbp;
    MPI_Comm_size(globCom, &nbp);
    int rank;
    MPI_Comm_rank(globCom, &rank);

    int global_size = params.discretization * params.discretization;
    
    if (rank == 0) {
        // Processus d'affichage (maître global)
        auto displayer = Displayer::init_instance(params.discretization, params.discretization);
        // Vecteurs pour stocker les cartes globales (feu et végétation)
        std::vector<unsigned char> global_fire_map(global_size);
        std::vector<unsigned char> global_veg_map(global_size);
        
        bool loop = true;
        while(loop) {
            int request = 1;
            // Demande de mise à jour des données au processus de calcul (rang 1)
            MPI_Send(&request, 1, MPI_INT, 1, 0, globCom);
            
            // Réception des cartes globales depuis le processus de calcul
            MPI_Recv(global_fire_map.data(), global_size, MPI_UNSIGNED_CHAR, 1, 0, globCom, MPI_STATUS_IGNORE);
            MPI_Recv(global_veg_map.data(), global_size, MPI_UNSIGNED_CHAR, 1, 1, globCom, MPI_STATUS_IGNORE);
            
            // Mise à jour de l'affichage via SDL
            displayer->update(global_veg_map, global_fire_map);
            
            // Gestion des événements SDL (fermeture de la fenêtre)
            SDL_Event event;
            if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
                int quit_signal = -1;
                MPI_Send(&quit_signal, 1, MPI_INT, 1, 0, globCom);
                loop = false;
            }
            std::this_thread::sleep_for(100ms);
        }
        SDL_Quit();
    }
    else {
        // Processus de calcul (simulation)
        Model simu(params.length, params.discretization, params.wind, params.start);
        bool continue_simulation = true;
        
        // Variables pour la mesure du temps
        double total_time = 0.0;
        int iter_count = 0;
        
        while (continue_simulation && simu.update()) {
            // Début de l'itération : mesure du temps
            double iter_start = MPI_Wtime();
            
            // Vérification non bloquante d'un éventuel signal d'arrêt venant du processus d'affichage
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
            
            // Fin de l'itération : calcul du temps écoulé
            double iter_end = MPI_Wtime();
            total_time += (iter_end - iter_start);
            iter_count++;
            
            // Pause optionnelle pour stabiliser l'affichage
            std::this_thread::sleep_for(100ms);
        }
        
        // Calcul et affichage du temps moyen par itération (calcul + communication)
        double average_time = (iter_count > 0) ? (total_time / iter_count) : 0.0;
        std::cout << "Temps moyen par itération (calcul + communication) : " 
                  << average_time << " secondes" << std::endl;
    }
    
    MPI_Finalize();
    return EXIT_SUCCESS;
}
