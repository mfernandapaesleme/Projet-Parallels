// simulation_parallèle.cpp
#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>
#include <mpi.h>
#include <omp.h>
#include <fstream>

#include "model.hpp"
#include "display.hpp"

using namespace std::string_literals;
using namespace std::chrono_literals;

// Structure pour stocker les paramètres de simulation
struct ParamsType {
    double length {1.0};             // Longueur du terrain (en km)
    unsigned discretization {300u};  // Nombre de cases par direction
    std::array<double,2> wind {0.0, 0.0};  // Vecteur vitesse du vent
    Model::LexicoIndices start {10u,10u};   // Position initiale du foyer (colonne, ligne)
    int num_threads {1};                    // Nombre de threads OpenMP
};

// Fonction récursive d'analyse des arguments de la ligne de commande
void analyze_arg(int nargs, char* args[], ParamsType& params) {
    if (nargs == 0) return;
    std::string key(args[0]);
    
    if (key == "-l"s || key.find("--longueur=") != std::string::npos) {
        if (nargs < 2) {
            std::cerr << "Manque une valeur pour la longueur du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        if (key == "-l"s) {
            params.length = std::stod(args[1]);
            analyze_arg(nargs-2, &args[2], params);
            return;
        }
        else {
            auto subkey = key.substr(key.find("=") + 1);
            params.length = std::stod(subkey);
            analyze_arg(nargs-1, &args[1], params);
            return;
        }
    }
    
    if (key == "-n"s || key.find("--number_of_cases=") != std::string::npos) {
        if (nargs < 2) {
            std::cerr << "Manque une valeur pour le nombre de cases par direction !" << std::endl;
            exit(EXIT_FAILURE);
        }
        if (key == "-n"s) {
            params.discretization = std::stoul(args[1]);
            analyze_arg(nargs-2, &args[2], params);
            return;
        }
        else {
            auto subkey = key.substr(key.find("=") + 1);
            params.discretization = std::stoul(subkey);
            analyze_arg(nargs-1, &args[1], params);
            return;
        }
    }
    
    if (key == "-w"s || key.find("--wind=") != std::string::npos) {
        if (nargs < 2) {
            std::cerr << "Manque une paire de valeurs pour la direction du vent !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values = (key == "-w"s) ? args[1] : key.substr(key.find("=") + 1);
        auto commaPos = values.find(",");
        if (commaPos == std::string::npos) {
            std::cerr << "Deux valeurs séparées par une virgule sont requises pour définir le vent" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.wind[0] = std::stod(values.substr(0, commaPos));
        params.wind[1] = std::stod(values.substr(commaPos + 1));
        if (key == "-w"s) {
            analyze_arg(nargs-2, &args[2], params);
        }
        return;
    }
    
    if (key == "-s"s || key.find("--start=") != std::string::npos) {
        if (nargs < 2) {
            std::cerr << "Manque une paire de valeurs pour la position initiale du foyer !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values = (key == "-s"s) ? args[1] : key.substr(key.find("=") + 1);
        auto commaPos = values.find(",");
        if (commaPos == std::string::npos) {
            std::cerr << "Deux valeurs séparées par une virgule sont requises pour définir la position initiale du foyer" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.start.column = std::stoul(values.substr(0, commaPos));
        params.start.row = std::stoul(values.substr(commaPos + 1));
        if (key == "-s"s) {
            analyze_arg(nargs-2, &args[2], params);
        }
        return;
    }
    
    // Option pour le nombre de threads
    if (key == "-t"s || key == "--threads"s) {
        if (nargs < 2) {
            std::cerr << "Manque une valeur pour le nombre de threads OpenMP !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.num_threads = std::stoi(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    
    // Si l'option n'est pas reconnue
    std::cerr << "Option inconnue : " << key << std::endl;
    exit(EXIT_FAILURE);
}

// Fonction de parsing des arguments
ParamsType parse_arguments(int nargs, char* args[]) {
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
        std::cerr << "[ERREUR FATALE] Le nombre de cellules par direction doit être positif et non nul !" << std::endl;
        flag = false;
    }
    if (params.start.row >= params.discretization || params.start.column >= params.discretization) {
        std::cerr << "[ERREUR FATALE] Mauvais indices pour la position initiale du foyer" << std::endl;
        flag = false;
    }
    if (params.num_threads <= 0) {
        std::cerr << "[ERREUR FATALE] Le nombre de threads doit être positif et non nul !" << std::endl;
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
              << params.start.column << ", " << params.start.row << std::endl
              << "\tNombre de threads OpenMP : " << params.num_threads << std::endl;
}

int main(int argc, char* argv[]) {
    // Traitement des arguments de la ligne de commande
    ParamsType params = parse_arguments(argc - 1, &argv[1]);
    display_params(params);
    if (!check_params(params)) {
        return EXIT_FAILURE;
    }

     // Initialisation de MPI avec support pour les threads multiples
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    
    if (provided < MPI_THREAD_MULTIPLE) {
        std::cerr << "MPI ne supporte pas le niveau de threading requis." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // Récupérer le nombre de processus et le rang
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Initialisation du nombre de threads OpenMP
    omp_set_num_threads(std::min(params.num_threads, omp_get_max_threads()));
    
    // Affichage des informations sur le matériel
    if (rank == 0) {
        std::cout << "Information sur l'hardware:" << std::endl
                  << "\tNombre maximum de threads OpenMP disponibles : " << omp_get_max_threads() << std::endl
                  << "\tNombre de processeurs physiques : " << omp_get_num_procs() << std::endl;
    }
    
    // Taille globale du terrain
    int global_size = params.discretization * params.discretization;
    
    if (rank == 0) {
        // Processus d'affichage (maître global)
        auto displayer = Displayer::init_instance(params.discretization, params.discretization);
        // Vecteurs pour stocker les cartes globales (feu et végétation)
        std::vector<unsigned char> global_fire_map(global_size, 0);
        std::vector<unsigned char> global_veg_map(global_size, 255);
        
        bool loop = true;
        while(loop) {
            int request = 1;
            // Demande de mise à jour des données au processus de calcul (rang 1)
            MPI_Send(&request, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);

            std::cout << "[Master] Aguardando atualização ou término..." << std::endl;
            MPI_Status status;
            int flag = 0;
            MPI_Iprobe(1, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
            std::cout << "[Master] Mensagem recebida do rank " << status.MPI_SOURCE << std::endl;

            if (status.MPI_TAG == 1 || status.MPI_TAG == 2 ) { // Mensagem de atualização
                std::cout << "[Master] Recebendo dados de atualização." << std::endl;
                MPI_Recv(global_fire_map.data(), global_size, MPI_UNSIGNED_CHAR, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(global_veg_map.data(), global_size, MPI_UNSIGNED_CHAR, 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::cout << "[Master] Mensagem de atualização recebida! " << std::endl;
                displayer->update(global_veg_map, global_fire_map);
            }
            else if (status.MPI_TAG == 3) { // Mensagem de término
                std::cout << "[Master] Recebendo sinal de término." << std::endl;
                MPI_Recv(NULL, 0, MPI_BYTE, 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::cout << "[Master] Mensagem de termino recebida! " << std::endl;
                loop = false; // Agora garantimos que o loop vai parar
            }
            
            // Gestion des événements SDL (fermeture de la fenêtre)
            SDL_Event event;
            if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
                int quit_signal = -1;
                MPI_Send(&quit_signal, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
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
        double compute_time = 0.0;  // Temps spécifique au calcul
        int iter_count = 0;
        
        std::cout << "Processus de calcul " << rank << " utilisant " << params.num_threads << " threads OpenMP." << std::endl;
        
        // Ouverture du fichier de log pour les processus de calcul
        std::string log_filename = "parallel_" + std::to_string(params.num_threads) + "_" + std::to_string(rank) + ".log";
        std::ofstream log_file(log_filename);
        log_file << "TimeStep Temps_Update Temps_Display Temps_Total\n";

        // Mesure du temps total de la simulation
        double start_total = MPI_Wtime();

        while (continue_simulation) {
            double iter_start = MPI_Wtime();

            // Mesure du temps d'Update
            double start_update = MPI_Wtime();
            // Avancement du feu déjà parallélisé avec OpenMP
            bool continue_simulation = simu.update();
            double end_update = MPI_Wtime();
            double elapsed_update = end_update - start_update;

            if (!continue_simulation) {
                std::cout << "[Calc " << rank << "] Enviando sinal de término para o mestre." << std::endl;
                MPI_Send(NULL, 0, MPI_BYTE, 0, 3, MPI_COMM_WORLD); // Tag=2 para término
            }
            
            // Vérification non bloquante d'un éventuel signal d'arrêt venant du processus d'affichage
            int flag = 0;
            MPI_Status status;
            MPI_Iprobe(0, 0, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                int signal;
                std::cout << "[Calc " << rank << "] Esperando receber algo" << std::endl;
                MPI_Recv(&signal, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::cout << "[Calc " << rank << "]  Mensagem algo recebida! " << std::endl;
                if (signal == -1) {
                    continue_simulation = false;
                    break;
                }
            }

            // Envoi des cartes de simulation (feu et végétation) au processus d'affichage
            std::cout << "[Calc " << rank << "] Enviando sinais de atualização para o mestre." << std::endl;
            MPI_Send(simu.fire_map().data(), global_size, MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD);
            MPI_Send(simu.vegetal_map().data(), global_size, MPI_UNSIGNED_CHAR, 0, 2, MPI_COMM_WORLD);
            
            // Fin de l'itération : calcul du temps écoulé
            double iter_end = MPI_Wtime();
            double elapsed_total = iter_end - start_total;
            total_time += (iter_end - iter_start);

            // Enregistrement des temps dans le fichier de log
            log_file << simu.time_step() << " " 
                     << elapsed_update << " " 
                     << elapsed_total << "\n";

            iter_count++;
            
            // Pause optionnelle pour stabiliser l'affichage
            std::this_thread::sleep_for(100ms);
        }
        
        // Calcul et affichage des statistiques de temps
        /* if (iter_count > 0) {
            double average_total_time = total_time / iter_count;
            double average_compute_time = compute_time / iter_count;
            double compute_percentage = (compute_time / total_time) * 100.0;
            
            std::cout << "=== Statistiques de performance ===" << std::endl;
            std::cout << "Nombre de threads OpenMP : " << params.num_threads << std::endl;
            std::cout << "Nombre d'itérations : " << iter_count << std::endl;
            std::cout << "Temps total moyen par itération : " << average_total_time << " secondes" << std::endl;
            std::cout << "Temps de calcul moyen par itération : " << average_compute_time << " secondes" << std::endl;
            std::cout << "Pourcentage du temps passé en calcul : " << compute_percentage << "%" << std::endl;
        } */
    }
    
    
    // Finalisation de MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    std::cout << "MPI finalizado" << std::endl;
    return EXIT_SUCCESS;
}