#include "model4.hpp"
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include <cstdlib>

// Função auxiliar para pseudo-aleatoriedade (simplificada)
namespace {
    double pseudo_random(std::size_t index, std::size_t time_step) {
        std::uint_fast32_t xi = std::uint_fast32_t(index*(time_step+1));
        std::uint_fast32_t r  = (48271 * xi) % 2147483647;
        return r / 2147483646.0;
    }
    double log_factor(std::uint8_t value) {
        return std::log(1. + value) / std::log(256);
    }
}

ModelMPI::ModelMPI(double t_length, unsigned dim, std::array<double,2> t_wind,
                   Model::LexicoIndices t_start, int rank, int size, double t_max_wind)
    : m_length(t_length), global_dim(dim), m_wind(t_wind), m_max_wind(t_max_wind),
      m_time_step(0), m_rank(rank), m_size(size)
{
    // Para uma grade quadrada, usamos "dim" para linhas e colunas
    // Dividimos as linhas entre os processos (deve ser divisível por m_size)
    local_dim = global_dim / m_size;
    local_grid_dim = local_dim + 2; // adiciona ghost cells (uma linha acima e outra abaixo)

    // Inicializa as grades locais
    vegetation.resize(local_grid_dim * global_dim, 255u);
    fire.resize(local_grid_dim * global_dim, 0u);

    // Inicializa o "fire front" local (como simplificação, usa o mapa de fogo)
    m_fire_front.clear();

    // Inicializa os parâmetros do modelo (valores simplificados)
    m_wind_speed = std::sqrt(m_wind[0]*m_wind[0] + m_wind[1]*m_wind[1]);
    constexpr double alpha0 = 4.52790762e-01;
    constexpr double alpha1 = 9.58264437e-04;
    constexpr double alpha2 = 3.61499382e-05;
    if (m_wind_speed < m_max_wind)
        p1 = alpha0 + alpha1 * m_wind_speed + alpha2 * (m_wind_speed * m_wind_speed);
    else 
        p1 = alpha0 + alpha1 * m_max_wind + alpha2 * (m_max_wind * m_max_wind);
    p2 = 0.3;
    
    if (m_wind[0] > 0) {
        alphaEastWest = std::abs(m_wind[0] / m_max_wind) + 1;
        alphaWestEast = 1 - std::abs(m_wind[0] / m_max_wind);
    } else {
        alphaWestEast = std::abs(m_wind[0] / m_max_wind) + 1;
        alphaEastWest = 1 - std::abs(m_wind[0] / m_max_wind);
    }
    if (m_wind[1] > 0) {
        alphaSouthNorth = std::abs(m_wind[1] / m_max_wind) + 1;
        alphaNorthSouth = 1 - std::abs(m_wind[1] / m_max_wind);
    } else {
        alphaNorthSouth = std::abs(m_wind[1] / m_max_wind) + 1;
        alphaSouthNorth = 1 - std::abs(m_wind[1] / m_max_wind);
    }

    // Inicializa o fogo no ponto inicial, se pertencer a este processo
    int global_fire_row = t_start.row;
    int global_fire_col = t_start.column;
    int owner = global_fire_row / local_dim;
    if (m_rank == owner) {
        // Converte para índice local (linha 0 é ghost)
        int local_fire_row = global_fire_row - owner * local_dim + 1;
        std::size_t index = local_fire_row * global_dim + global_fire_col;
        fire[index] = 255u;
        m_fire_front[index] = 255u;
    }
}

// Troca as ghost cells para o vetor grid (vegetation ou fire)
void ModelMPI::exchangeGhostRows(std::vector<std::uint8_t>& grid) {
    MPI_Status status;
    // Troca com o processo acima
    if (m_rank > 0) {
        MPI_Sendrecv(&grid[global_dim], global_dim, MPI_UNSIGNED_CHAR, m_rank - 1, 0,
                     &grid[0], global_dim, MPI_UNSIGNED_CHAR, m_rank - 1, 0,
                     MPI_COMM_WORLD, &status);
    }
    // Troca com o processo abaixo
    if (m_rank < m_size - 1) {
        MPI_Sendrecv(&grid[local_dim * global_dim], global_dim, MPI_UNSIGNED_CHAR, m_rank + 1, 0,
                     &grid[(local_dim + 1) * global_dim], global_dim, MPI_UNSIGNED_CHAR, m_rank + 1, 0,
                     MPI_COMM_WORLD, &status);
    }
}

// Atualiza a simulação para o subdomínio deste processo
bool ModelMPI::update() {
    exchangeGhostRows(vegetation);
    exchangeGhostRows(fire);
    updateLocalSimulation();
    m_time_step++;
    return !m_fire_front.empty(); // retorna se ainda há fogo ativo localmente
}

// Atualiza a simulação local (exclui ghost cells) com regra simplificada de propagação
void ModelMPI::updateLocalSimulation() {
    std::vector<std::uint8_t> new_fire = fire; // Cópia para atualização síncrona
    m_fire_front.clear();

    // Percorre as linhas reais (de 1 a local_dim)
    for (unsigned i = 1; i <= local_dim; i++) {
        for (unsigned j = 0; j < global_dim; j++) {
            std::size_t idx = i * global_dim + j;
            if (fire[idx] == 255) {
                // Propaga para os vizinhos (Norte, Sul, Oeste, Leste)
                if (i - 1 >= 0) {
                    std::size_t nidx = (i - 1) * global_dim + j;
                    if (vegetation[nidx] > 0 && fire[nidx] == 0)
                        new_fire[nidx] = 255;
                }
                if (i + 1 < local_grid_dim) {
                    std::size_t nidx = (i + 1) * global_dim + j;
                    if (vegetation[nidx] > 0 && fire[nidx] == 0)
                        new_fire[nidx] = 255;
                }
                if (j > 0) {
                    std::size_t nidx = i * global_dim + (j - 1);
                    if (vegetation[nidx] > 0 && fire[nidx] == 0)
                        new_fire[nidx] = 255;
                }
                if (j + 1 < global_dim) {
                    std::size_t nidx = i * global_dim + (j + 1);
                    if (vegetation[nidx] > 0 && fire[nidx] == 0)
                        new_fire[nidx] = 255;
                }
                // Simula a diminuição do fogo
                new_fire[idx] = 128;
            } else if (fire[idx] != 0) {
                new_fire[idx] = fire[idx] >> 1;
            }
            if (new_fire[idx] > 0)
                m_fire_front[idx] = new_fire[idx];
        }
    }
    fire = new_fire;

    // Diminui a vegetação onde há fogo
    for (auto& f : m_fire_front) {
        std::size_t idx = f.first;
        if (vegetation[idx] > 0)
            vegetation[idx] -= 1;
    }
}

bool ModelMPI::checkLocalFire() const {
    for (unsigned i = 1; i <= local_dim; i++) {
        for (unsigned j = 0; j < global_dim; j++) {
            if (fire[i * global_dim + j] > 0)
                return true;
        }
    }
    return false;
}

bool ModelMPI::localFireActive() const {
    return checkLocalFire();
}

void ModelMPI::getLocalData(std::vector<std::uint8_t>& localVegetation,
                            std::vector<std::uint8_t>& localFire) const {
    localVegetation.resize(local_dim * global_dim);
    localFire.resize(local_dim * global_dim);
    for (unsigned i = 1; i <= local_dim; i++) {
        for (unsigned j = 0; j < global_dim; j++) {
            localVegetation[(i - 1) * global_dim + j] = vegetation[i * global_dim + j];
            localFire[(i - 1) * global_dim + j] = fire[i * global_dim + j];
        }
    }
}
