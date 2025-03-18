#include <stdexcept>
#include <cmath>
#include <iostream>
#include <vector>
#include <omp.h>
#include "model.hpp"


namespace
{
    double pseudo_random( std::size_t index, std::size_t time_step )
    {
        std::uint_fast32_t xi = std::uint_fast32_t(index*(time_step+1));
        std::uint_fast32_t r  = (48271*xi)%2147483647;
        return r/2147483646.;
    }

    double log_factor( std::uint8_t value )
    {
        return std::log(1.+value)/std::log(256);
    }
}

Model::Model( double t_length, unsigned t_discretization, std::array<double,2> t_wind,
              LexicoIndices t_start_fire_position, double t_max_wind )
    :   m_length(t_length),
        m_distance(-1),
        m_geometry(t_discretization),
        m_wind(t_wind),
        m_wind_speed(std::sqrt(t_wind[0]*t_wind[0] + t_wind[1]*t_wind[1])),
        m_max_wind(t_max_wind),
        m_vegetation_map(t_discretization*t_discretization, 255u),
        m_fire_map(t_discretization*t_discretization, 0u),
        m_time_step(0)
{
    if (t_discretization == 0)
    {
        throw std::range_error("Le nombre de cases par direction doit être plus grand que zéro.");
    }
    m_distance = m_length/double(m_geometry);
    auto index = get_index_from_lexicographic_indices(t_start_fire_position);
    m_fire_map[index] = 255u;
    m_fire_front[index] = 255u;

    constexpr double alpha0 = 4.52790762e-01;
    constexpr double alpha1 = 9.58264437e-04;
    constexpr double alpha2 = 3.61499382e-05;

    if (m_wind_speed < t_max_wind)
        p1 = alpha0 + alpha1*m_wind_speed + alpha2*(m_wind_speed*m_wind_speed);
    else 
        p1 = alpha0 + alpha1*t_max_wind + alpha2*(t_max_wind*t_max_wind);
    p2 = 0.3;

    if (m_wind[0] > 0)
    {
        alphaEastWest = std::abs(m_wind[0]/t_max_wind)+1;
        alphaWestEast = 1.-std::abs(m_wind[0]/t_max_wind);    
    }
    else
    {
        alphaWestEast = std::abs(m_wind[0]/t_max_wind)+1;
        alphaEastWest = 1. - std::abs(m_wind[0]/t_max_wind);
    }

    if (m_wind[1] > 0)
    {
        alphaSouthNorth = std::abs(m_wind[1]/t_max_wind) + 1;
        alphaNorthSouth = 1. - std::abs(m_wind[1]/t_max_wind);
    }
    else
    {
        alphaNorthSouth = std::abs(m_wind[1]/t_max_wind) + 1;
        alphaSouthNorth = 1. - std::abs(m_wind[1]/t_max_wind);
    }
}

bool Model::update()
{
    // Extrair todas as chaves do m_fire_front para um vetor
    std::vector<std::size_t> fire_front_keys;
    fire_front_keys.reserve(m_fire_front.size());
    for (const auto& f : m_fire_front) {
        fire_front_keys.push_back(f.first);
    }

    // Inicializar vetores de vetores para armazenar contribuições de cada thread
    std::vector<std::vector<std::size_t>> new_fire_cells_threads;
    std::vector<std::vector<std::pair<std::size_t, std::uint8_t>>> updated_cells_threads;

    int num_threads;
    #pragma omp parallel
    {
        #pragma omp single
        {
            num_threads = omp_get_num_threads();
            new_fire_cells_threads.resize(num_threads);
            updated_cells_threads.resize(num_threads);
        }

        int thread_id = omp_get_thread_num();
        // Garantir que thread_id está dentro do range
        if (thread_id < new_fire_cells_threads.size()) {
            std::vector<std::size_t>& local_new_fire = new_fire_cells_threads[thread_id];
            std::vector<std::pair<std::size_t, std::uint8_t>>& local_updates = updated_cells_threads[thread_id];

            #pragma omp for schedule(static)
            for (size_t i = 0; i < fire_front_keys.size(); i++) {
                std::size_t key = fire_front_keys[i];
                // Verificar se key é válido
                if (key >= m_fire_map.size() || key >= m_vegetation_map.size()) {
                    std::cerr << "Erro: key " << key << " está fora dos limites.\n";
                    continue;
                }

                std::uint8_t value = m_fire_front[key];

                // Recuperação da coordenada lexicográfica da célula em fogo
                LexicoIndices coord = get_lexicographic_from_index(key);

                // Validar coordenadas
                if (coord.row >= m_geometry || coord.column >= m_geometry) {
                    std::cerr << "Erro: coordenada (row: " << coord.row << ", column: " << coord.column << ") inválida.\n";
                    continue;
                }

                double power = log_factor(value);

                // Verificar células vizinhas (sul, norte, leste, oeste)

                // Sul
                if (coord.row < m_geometry-1) {
                    std::size_t south_key = key + m_geometry;
                    if (south_key >= m_fire_map.size()) {
                        std::cerr << "Erro: south_key " << south_key << " está fora dos limites.\n";
                        continue;
                    }
                    double tirage = pseudo_random(key+m_time_step, m_time_step);
                    double green_power = m_vegetation_map[south_key];
                    double correction = power * log_factor(green_power);
                    if (tirage < alphaSouthNorth*p1*correction && m_fire_map[south_key] == 0) {
                        local_new_fire.push_back(south_key);
                    }
                }

                // Norte
                if (coord.row > 0) {
                    std::size_t north_key = key - m_geometry;
                    if (north_key >= m_fire_map.size()) {
                        std::cerr << "Erro: north_key " << north_key << " está fora dos limites.\n";
                        continue;
                    }
                    double tirage = pseudo_random(key*13427+m_time_step, m_time_step);
                    double green_power = m_vegetation_map[north_key];
                    double correction = power * log_factor(green_power);
                    if (tirage < alphaNorthSouth*p1*correction && m_fire_map[north_key] == 0) {
                        local_new_fire.push_back(north_key);
                    }
                }

                // Leste
                if (coord.column < m_geometry-1) {
                    std::size_t east_key = key + 1;
                    if (east_key >= m_fire_map.size()) {
                        std::cerr << "Erro: east_key " << east_key << " está fora dos limites.\n";
                        continue;
                    }
                    double tirage = pseudo_random(key*13427*13427+m_time_step, m_time_step);
                    double green_power = m_vegetation_map[east_key];
                    double correction = power * log_factor(green_power);
                    if (tirage < alphaEastWest*p1*correction && m_fire_map[east_key] == 0) {
                        local_new_fire.push_back(east_key);
                    }
                }

                // Oeste
                if (coord.column > 0) {
                    std::size_t west_key = key - 1;
                    if (west_key >= m_fire_map.size()) {
                        std::cerr << "Erro: west_key " << west_key << " está fora dos limites.\n";
                        continue;
                    }
                    double tirage = pseudo_random(key*13427*13427*13427+m_time_step, m_time_step);
                    double green_power = m_vegetation_map[west_key];
                    double correction = power * log_factor(green_power);
                    if (tirage < alphaWestEast*p1*correction && m_fire_map[west_key] == 0) {
                        local_new_fire.push_back(west_key);
                    }
                }

                // Atualizando o estado atual do fogo
                std::uint8_t new_value = value;
                if (value == 255) {
                    double tirage = pseudo_random(key * 52513 + m_time_step, m_time_step);
                    if (tirage < p2) {
                        new_value = value >> 1; // Divide por 2
                        local_updates.emplace_back(key, new_value);
                    }
                } else {
                    new_value = value >> 1; // Divide por 2
                    local_updates.emplace_back(key, new_value);
                }
            }
        }
    } // Fim da região paralela

    // Combinar os resultados de todas as threads fora da região paralela
    std::vector<std::size_t> new_fire_cells;
    std::vector<std::pair<std::size_t, std::uint8_t>> updated_cells;

    // Estimar o tamanho total e reservar espaço para evitar realocações
    size_t total_new_fire = 0;
    size_t total_updated = 0;
    for (int t = 0; t < num_threads; ++t) {
        total_new_fire += new_fire_cells_threads[t].size();
        total_updated += updated_cells_threads[t].size();
    }
    new_fire_cells.reserve(total_new_fire);
    updated_cells.reserve(total_updated);

    for (int t = 0; t < num_threads; ++t) {
        new_fire_cells.insert(new_fire_cells.end(), new_fire_cells_threads[t].begin(), new_fire_cells_threads[t].end());
        updated_cells.insert(updated_cells.end(), updated_cells_threads[t].begin(), updated_cells_threads[t].end());
    }

    // **Opcional**: Remover duplicatas, se necessário
    // Exemplo: se múltiplas threads adicionaram o mesmo `south_key`
    // std::sort(new_fire_cells.begin(), new_fire_cells.end());
    // new_fire_cells.erase(std::unique(new_fire_cells.begin(), new_fire_cells.end()), new_fire_cells.end());

    // Aplicar as mudanças nas células globais
    for (const auto& key : new_fire_cells) {
        if (key < m_fire_map.size()) {
            m_fire_map[key] = 255;
            m_fire_front[key] = 255;
        } else {
            std::cerr << "Erro: new_fire_cells key " << key << " está fora dos limites do m_fire_map.\n";
        }
    }

    for (const auto& update : updated_cells) {
        std::size_t key = update.first;
        std::uint8_t new_value = update.second;

        if (key < m_fire_map.size()) {
            m_fire_map[key] = new_value;
            if (new_value == 0) {
                m_fire_front.erase(key);
            } else {
                m_fire_front[key] = new_value;
            }
        } else {
            std::cerr << "Erro: updated_cells key " << key << " está fora dos limites do m_fire_map.\n";
        }
    }

    // Atualização da vegetação - mantida paralelizada
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < fire_front_keys.size(); i++) {
        std::size_t key = fire_front_keys[i];
        if (key < m_vegetation_map.size() && m_vegetation_map[key] > 0) {
            #pragma omp atomic
            m_vegetation_map[key] -= 1;
        } else {
            // Opcional: Log se a chave está fora dos limites ou se a vegetação já está 0
            // std::cerr << "Aviso: key " << key << " fora dos limites ou vegetação já 0.\n";
        }
    }

    m_time_step += 1;
    return !m_fire_front.empty();
}

std::size_t   
Model::get_index_from_lexicographic_indices( LexicoIndices t_lexico_indices  ) const
{
    return t_lexico_indices.row*this->geometry() + t_lexico_indices.column;
}

auto 
Model::get_lexicographic_from_index( std::size_t t_global_index ) const -> LexicoIndices
{
    LexicoIndices ind_coords;
    ind_coords.row    = t_global_index/this->geometry();
    ind_coords.column = t_global_index%this->geometry();
    return ind_coords;
}