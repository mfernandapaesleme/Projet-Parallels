#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>
#include <omp.h>
#include <fstream>

// Inclua o cabeçalho da TBB para medir tempo
#include <tbb/tick_count.h>

#include "model.hpp"
#include "display.hpp"

using namespace std::string_literals;
using namespace std::chrono_literals;

struct ParamsType
{
    double length{1.};
    unsigned discretization{20u};
    std::array<double,2> wind{0.,0.};
    Model::LexicoIndices start{10u,10u};
    int num_threads{1};
};

void analyze_arg( int nargs, char* args[], ParamsType& params )
{
    if (nargs ==0) return;
    std::string key(args[0]);
    if (key == "-l"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour la longueur du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.length = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    auto pos = key.find("--longueur=");
    if (pos < key.size())
    {
        auto subkey = std::string(key,pos+11);
        params.length = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-n"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour le nombre de cases par direction pour la discrétisation du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.discretization = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--number_of_cases=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+18);
        params.discretization = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-w"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs pour la direction du vent !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values = std::string(args[1]);
        params.wind[0] = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--wind=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+7);
        params.wind[0] = std::stoul(subkey);
        auto pos2 = subkey.find(",");
        if (pos2 == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos2+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-s"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs para a position du foyer initial !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values = std::string(args[1]);
        params.start.column = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la position du foyer initial" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--start=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+8);
        params.start.column = std::stoul(subkey);
        auto pos2 = subkey.find(",");
        if (pos2 == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par uma virgule para definir a posição do fogo" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos2+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    // Nova opção para número de threads
    if (key == "-t"s || key == "--threads"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour le nombre de threads OpenMP !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.num_threads = std::stoi(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
}

ParamsType parse_arguments( int nargs, char* args[] )
{
    if (nargs == 0) return {};
    if ( (std::string(args[0]) == "--help"s) || (std::string(args[0]) == "-h") )
    {
        std::cout << 
R"RAW(Usage : simulation [option(s)]
  Lance la simulation d'incendie en prenant en compte les [option(s)].
  Les options sont :
    -l, --longueur=LONGUEUR     Définit la taille LONGUEUR (réel en km) du carré représentant la carte de la végétation.
    -n, --number_of_cases=N     Nombre n de cases par direction pour la discrétisation
    -w, --wind=VX,VY            Définit le vecteur vitesse du vent (pas de vent par défaut).
    -s, --start=COL,ROW         Définit les indices I,J de la case où commence l'incendie (milieu de la carte par défaut)
    -t, --threads=N             Définit le nombre de threads OpenMP à utiliser (défaut: 1)
    )RAW";
        exit(EXIT_SUCCESS);
    }
    ParamsType params;
    analyze_arg(nargs, args, params);
    return params;
}

bool check_params(ParamsType& params)
{
    bool flag = true;
    if (params.length <= 0)
    {
        std::cerr << "[ERREUR FATALE] La longueur du terrain doit être positive et non nulle !" << std::endl;
        flag = false;
    }

    if (params.discretization <= 0)
    {
        std::cerr << "[ERREUR FATALE] Le nombre de cellules par direction doit être positive et non nulle !" << std::endl;
        flag = false;
    }

    if ( (params.start.row >= params.discretization) || (params.start.column >= params.discretization) )
    {
        std::cerr << "[ERREUR FATALE] Mauvais indices pour la position initiale du foyer" << std::endl;
        flag = false;
    }

    if (params.num_threads <= 0)
    {
        std::cerr << "[ERREUR FATALE] Le nombre de threads doit être positif et non nul !" << std::endl;
        flag = false;
    }
    
    return flag;
}

void display_params(ParamsType const& params)
{
    std::cout << "Parametres définis pour la simulation : \n"
              << "\tTaille du terrain : " << params.length << std::endl 
              << "\tNombre de cellules par direction : " << params.discretization << std::endl 
              << "\tVecteur vitesse : [" << params.wind[0] << ", " << params.wind[1] << "]" << std::endl
              << "\tPosition initiale du foyer (col, ligne) : " << params.start.column << ", " << params.start.row << std::endl
              << "\tNombre de threads OpenMP : " << params.num_threads << std::endl;
}

int main( int nargs, char* args[] )
{
    // 1) Lê e valida parâmetros
    auto params = parse_arguments(nargs-1, &args[1]);
    display_params(params);
    if (!check_params(params)) return EXIT_FAILURE;

    // 2) Define o número de threads para OpenMP (independente da TBB)
    omp_set_num_threads(std::min(params.num_threads, omp_get_max_threads()));
    
    // Informações sobre hardware
    std::cout << "Information sur l'hardware:" << std::endl;
    std::cout << "\tNúmero máximo de threads OpenMP disponibles: " << omp_get_max_threads() << std::endl;
    std::cout << "\tNúmero de processadores physiques: " << omp_get_num_procs() << std::endl;

    // 3) Cria o Displayer e o Model
    auto displayer = Displayer::init_instance( params.discretization, params.discretization );
    auto simu = Model( params.length, 
                       params.discretization, 
                       params.wind,
                       params.start );

    SDL_Event event;

    // 4) Abre arquivo de saída para registrar os tempos
    std::ofstream file("temps.txt");
    file << "TimeStep Temps_Avancement Temps_Affichage Temps_Simulation\n"; 

    // 5) Marcador de tempo inicial (TBB)
    tbb::tick_count start_total = tbb::tick_count::now();

    // 6) Loop principal de simulação
    while (true)
    {
        // Processa eventos (fechar janela, etc.)
        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_QUIT)
                return EXIT_SUCCESS;
        }

        // ---------------------------
        // a) Tempo de Avanço (update)
        // ---------------------------
        tbb::tick_count start_update = tbb::tick_count::now();
        bool isrunning = simu.update();  
        double Temps_Avancement = (tbb::tick_count::now() - start_update).seconds();

        if(!isrunning) {
            // Se o fogo acabou, sai do loop
            break;
        }

        // Exemplo de exibir o time_step a cada 32 iterações
        if ((simu.time_step() & 31) == 0) {
            std::cout << "Time step " << simu.time_step() << "\n===============" << std::endl;
        }

        // ------------------------------
        // b) Tempo de Exibição (display)
        // ------------------------------
        tbb::tick_count start_display = tbb::tick_count::now();
        displayer->update(simu.vegetal_map(), simu.fire_map());
        double Temps_Affichage = (tbb::tick_count::now() - start_display).seconds();

        // (opcional) Pausa de 50ms para não “atropelar” a CPU
        std::this_thread::sleep_for(50ms);

        // ------------------------------
        // c) Tempo total do Time Step
        // ------------------------------
        double Temps_Simulation = (tbb::tick_count::now() - start_total).seconds();

        // 7) Registra os tempos no arquivo
        file << simu.time_step()        << " "
             << Temps_Avancement       << " "
             << Temps_Affichage        << " "
             << Temps_Simulation       << "\n";
    }

    file.close();

    // 8) Tempo total da simulação
    double elapsed_total = (tbb::tick_count::now() - start_total).seconds();
    std::cout << "Temps total de la simulation: " << elapsed_total << " secondes\n";

    return EXIT_SUCCESS;
}
