#ifndef MODEL_MPI_HPP
#define MODEL_MPI_HPP

#include <cstdint>
#include <array>
#include <vector>
#include <unordered_map>
#include "model.hpp"  // Reutiliza a definição de ParamsType e de LexicoIndices

/**
 * A classe ModelMPI encapsula toda a lógica da simulação paralela.
 * Supõe que a grade é quadrada, de dimensão "dim" x "dim".
 */
class ModelMPI {
public:
    struct LexicoIndices {
        std::size_t row, column;
        LexicoIndices(std::size_t t_column = 0, std::size_t t_row = 0)
            : row(t_row), column(t_column) {}
    };

    ModelMPI(double t_length, unsigned dim, std::array<double,2> t_wind,
             Model::LexicoIndices t_start, int rank, int size, double t_max_wind = 60.);

    // Atualiza a simulação: troca ghost cells e aplica as regras de propagação
    bool update();

    // Retorna os dados locais (sem ghost cells) para visualização ou gathering
    void getLocalData(std::vector<std::uint8_t>& localVegetation,
                      std::vector<std::uint8_t>& localFire) const;

    // Verifica se há fogo ativo no subdomínio local (linhas reais)
    bool localFireActive() const;

    // Retorna a dimensão local (número de linhas reais) e a dimensão global
    unsigned getLocalDim() const { return local_dim; }
    unsigned getGlobalDim() const { return global_dim; }
    std::size_t time_step() const { return m_time_step; }

private:
    // Funções internas
    void exchangeGhostRows(std::vector<std::uint8_t>& grid);
    void updateLocalSimulation();
    bool checkLocalFire() const;

    // Parâmetros do domínio global
    double m_length;     // Tamanho do terreno global (em km)
    unsigned global_dim; // Dimensão global da grade (global_dim x global_dim)
    
    // Parâmetros do subdomínio local
    unsigned local_dim;       // Número de linhas reais deste processo
    unsigned local_grid_dim;  // local_dim + 2 (inclui ghost cells superior e inferior)
    
    // Dados locais (grade de vegetação e de fogo)
    std::vector<std::uint8_t> vegetation; // Inicializado com 255
    std::vector<std::uint8_t> fire;         // Inicializado com 0
    
    // Parâmetros do modelo (ex.: vento)
    std::array<double,2> m_wind;
    double m_wind_speed;
    double m_max_wind;
    double p1, p2;
    double alphaEastWest, alphaWestEast, alphaSouthNorth, alphaNorthSouth;

    // Gerenciamento do “fire front” (simplificado)
    std::unordered_map<std::size_t, std::uint8_t> m_fire_front;

    // Controle do tempo da simulação
    std::size_t m_time_step;

    // Informações MPI
    int m_rank;
    int m_size;
};

#endif // MODEL_MPI_HPP
