#include <mpi.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <thread>
#include <chrono>
#include <fstream>
#include "model4.hpp"    // Reutilizamos as funções de parse e check dos parâmetros
#include "display.hpp"  // Para visualização (SDL)

using namespace std::string_literals;
using namespace std::chrono_literals;

struct ParamsType
{
    double length{1.};
    unsigned discretization{300u};
    std::array<double,2> wind{0.,0.};
    Model::LexicoIndices start{10u,10u};
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
        std::string values =std::string(args[1]);
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
        auto pos = subkey.find(",");
        if (pos == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-s"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs pour la position du foyer initial !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values =std::string(args[1]);
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
        auto pos = subkey.find(",");
        if (pos == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params);
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
    
    return flag;
}

void display_params(ParamsType const& params)
{
    std::cout << "Parametres définis pour la simulation : \n"
              << "\tTaille du terrain : " << params.length << std::endl 
              << "\tNombre de cellules par direction : " << params.discretization << std::endl 
              << "\tVecteur vitesse : [" << params.wind[0] << ", " << params.wind[1] << "]" << std::endl
              << "\tPosition initiale du foyer (col, ligne) : " << params.start.column << ", " << params.start.row << std::endl;
}

// Função para trocar as linhas de fronteira (ghost cells) entre processos
// 'grid' representa o vetor (em 1D) que armazena a grade (vegetation ou fire)
// cols: número de colunas (tamanho horizontal da grade)
// local_rows: número de linhas “reais” (sem as linhas de ghost) do subdomínio
// Cada processo aloca (local_rows + 2) linhas: uma de ghost no topo e outra no final.
void exchangeGhostRows(std::vector<std::uint8_t>& grid, int cols, int local_rows, int rank, int size) {
    MPI_Status status;
    // Troca com o processo acima (se existir)
    if (rank > 0) {
        // Envia a primeira linha real (linha 1) para o processo acima e recebe a ghost row na linha 0
        MPI_Sendrecv(&grid[cols], cols, MPI_UNSIGNED_CHAR, rank - 1, 0,
                     &grid[0], cols, MPI_UNSIGNED_CHAR, rank - 1, 0,
                     MPI_COMM_WORLD, &status);
    }
    // Troca com o processo abaixo (se existir)
    if (rank < size - 1) {
        // Envia a última linha real (linha local_rows) para o processo abaixo e recebe a ghost row na linha local_rows+1
        MPI_Sendrecv(&grid[local_rows * cols], cols, MPI_UNSIGNED_CHAR, rank + 1, 0,
                     &grid[(local_rows + 1) * cols], cols, MPI_UNSIGNED_CHAR, rank + 1, 0,
                     MPI_COMM_WORLD, &status);
    }
}

// Uma versão simplificada do update da simulação para o subdomínio local.
// Essa função percorre as células reais (excluindo as linhas de ghost) e aplica uma regra simples de propagação do fogo.
void update_local_simulation(std::vector<std::uint8_t>& vegetation,
                             std::vector<std::uint8_t>& fire,
                             int cols, int local_rows) {
    std::vector<std::uint8_t> new_fire = fire; // Cópia para atualizar de forma síncrona

    // Percorre as linhas reais: de 1 até local_rows (a linha 0 e a linha local_rows+1 são ghost)
    for (int i = 1; i <= local_rows; i++) {
        for (int j = 0; j < cols; j++) {
            int idx = i * cols + j;
            if (fire[idx] == 255) { // Célula em chamas
                // Propaga para vizinhos (Norte, Sul, Leste e Oeste)
                // Norte
                if (i - 1 >= 0) {
                    int nidx = (i - 1) * cols + j;
                    if (vegetation[nidx] > 0 && fire[nidx] == 0)
                        new_fire[nidx] = 255;
                }
                // Sul
                if (i + 1 < local_rows + 2) {
                    int nidx = (i + 1) * cols + j;
                    if (vegetation[nidx] > 0 && fire[nidx] == 0)
                        new_fire[nidx] = 255;
                }
                // Oeste
                if (j - 1 >= 0) {
                    int nidx = i * cols + (j - 1);
                    if (vegetation[nidx] > 0 && fire[nidx] == 0)
                        new_fire[nidx] = 255;
                }
                // Leste
                if (j + 1 < cols) {
                    int nidx = i * cols + (j + 1);
                    if (vegetation[nidx] > 0 && fire[nidx] == 0)
                        new_fire[nidx] = 255;
                }
                // Simula a diminuição do fogo (aqui, simplificamos para um valor intermediário)
                new_fire[idx] = 128;
            }
            else if (fire[idx] != 0) {
                // Reduz o fogo (decai)
                new_fire[idx] = fire[idx] >> 1;
            }
        }
    }
    fire = new_fire;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Ouverture du fichier pour enregistrer les temps (seulement sur le processus 0)
    std::ofstream file;
    if (rank == 0) {
        file.open("temps.txt");
        file << "TimeStep Temps_Avancement Temps_Affichage Temps_Simulation\n";
    }
    
    // Traitement des arguments
    auto params = parse_arguments(argc - 1, &argv[1]);
    if (!check_params(params)) {
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    
    if (rank == 0) {
        display_params(params);
    }
    
    int global_dim = params.discretization;
    if (global_dim % size != 0) {
        if (rank == 0)
            std::cerr << "O número de processos deve dividir exatamente a dimensão (" 
                      << global_dim << ")." << std::endl;
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    
    // Création de l'instance du modèle MPI (ModelMPI)
    ModelMPI model(params.length, global_dim, params.wind, params.start, rank, size);
    
    std::shared_ptr<Displayer> displayer;
    if (rank == 0)
        displayer = Displayer::init_instance(global_dim, global_dim);
    
    bool global_continue = true;
    int time_step = 0;
    
    // Début de la mesure du temps total avec MPI_Wtime()
    double start_total = MPI_Wtime();
    
    while (global_continue) {
        double start_update = 0.0, end_update = 0.0, start_display = 0.0, end_display = 0.0;
        if (rank == 0)
            start_update = MPI_Wtime();
        
        model.update();
        
        if (rank == 0)
            end_update = MPI_Wtime();
        
        bool local_active = model.localFireActive();
        MPI_Allreduce(&local_active, &global_continue, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
        
        if (rank == 0) {
            std::vector<std::uint8_t> full_fire(global_dim * global_dim, 0);
            std::vector<std::uint8_t> full_vegetation(global_dim * global_dim, 0);
            
            std::vector<std::uint8_t> local_fire, local_vegetation;
            model.getLocalData(local_vegetation, local_fire);
            std::copy(local_fire.begin(), local_fire.end(), full_fire.begin());
            std::copy(local_vegetation.begin(), local_vegetation.end(), full_vegetation.begin());
            
            int local_size = model.getLocalDim() * model.getGlobalDim();
            MPI_Gather(MPI_IN_PLACE, local_size, MPI_UNSIGNED_CHAR,
                       full_fire.data(), local_size, MPI_UNSIGNED_CHAR,
                       0, MPI_COMM_WORLD);
            MPI_Gather(MPI_IN_PLACE, local_size, MPI_UNSIGNED_CHAR,
                       full_vegetation.data(), local_size, MPI_UNSIGNED_CHAR,
                       0, MPI_COMM_WORLD);
            
            start_display = MPI_Wtime();
            displayer->update(full_vegetation, full_fire);
            end_display = MPI_Wtime();
            
            double elapsed_update = end_update - start_update;
            double elapsed_display = end_display - start_display;
            double elapsed_total = MPI_Wtime() - start_total;
            
            file << time_step << " " 
                 << elapsed_update << " " 
                 << elapsed_display << " " 
                 << elapsed_total << "\n";
            
            std::cout << "Time step " << time_step << "\n===============" << std::endl;
        } else {
            int local_size = model.getLocalDim() * model.getGlobalDim();
            std::vector<std::uint8_t> local_fire, local_vegetation;
            model.getLocalData(local_vegetation, local_fire);
            MPI_Gather(local_fire.data(), local_size, MPI_UNSIGNED_CHAR,
                       nullptr, local_size, MPI_UNSIGNED_CHAR,
                       0, MPI_COMM_WORLD);
            MPI_Gather(local_vegetation.data(), local_size, MPI_UNSIGNED_CHAR,
                       nullptr, local_size, MPI_UNSIGNED_CHAR,
                       0, MPI_COMM_WORLD);
        }
        
        time_step++;
        std::this_thread::sleep_for(100ms);
    }
    
    if (rank == 0) {
        file.close();
    }
    
    MPI_Finalize();
    return EXIT_SUCCESS;
}