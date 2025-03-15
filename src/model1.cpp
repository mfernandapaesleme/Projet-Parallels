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

    // Criar vetores temporários para armazenar as novas células em chamas
    // e as células que se alteraram
    std::vector<std::size_t> new_fire_cells;
    std::vector<std::pair<std::size_t, std::uint8_t>> updated_cells;
    
    // Reservar espaço para evitar realocações
    new_fire_cells.reserve(fire_front_keys.size() * 4); // No máximo 4 vizinhos por célula
    updated_cells.reserve(fire_front_keys.size());

    // Processar cada célula em paralelo, coletando mudanças em vez de aplicá-las
    #pragma omp parallel
    {
        // Vetores locais para cada thread
        std::vector<std::size_t> local_new_fire;
        std::vector<std::pair<std::size_t, std::uint8_t>> local_updates;

        // Pré-alocação para evitar realocações dentro do loop
        local_new_fire.reserve(fire_front_keys.size());
        local_updates.reserve(fire_front_keys.size());

        #pragma omp for schedule(dynamic, 64)
        for (size_t i = 0; i < fire_front_keys.size(); i++) {
            std::size_t key = fire_front_keys[i];
            std::uint8_t value = m_fire_front[key];
            
            // Recuperação da coordenada lexicográfica da case em fogo
            LexicoIndices coord = get_lexicographic_from_index(key);
            // E da puissance do foyer
            double power = log_factor(value);

            // Verificar se a célula ao sul pode ser contaminada
            if (coord.row < m_geometry-1) {
                std::size_t south_key = key + m_geometry;
                double tirage = pseudo_random(key+m_time_step, m_time_step);
                double green_power = m_vegetation_map[south_key];
                double correction = power*log_factor(green_power);
                if (tirage < alphaSouthNorth*p1*correction && m_fire_map[south_key] == 0) {
                    local_new_fire.push_back(south_key);
                }
            }

            // Verificar se a célula ao norte pode ser contaminada
            if (coord.row > 0) {
                std::size_t north_key = key - m_geometry;
                double tirage = pseudo_random(key*13427+m_time_step, m_time_step);
                double green_power = m_vegetation_map[north_key];
                double correction = power*log_factor(green_power);
                if (tirage < alphaNorthSouth*p1*correction && m_fire_map[north_key] == 0) {
                    local_new_fire.push_back(north_key);
                }
            }

            // Verificar se a célula a leste pode ser contaminada
            if (coord.column < m_geometry-1) {
                std::size_t east_key = key + 1;
                double tirage = pseudo_random(key*13427*13427+m_time_step, m_time_step);
                double green_power = m_vegetation_map[east_key];
                double correction = power*log_factor(green_power);
                if (tirage < alphaEastWest*p1*correction && m_fire_map[east_key] == 0) {
                    local_new_fire.push_back(east_key);
                }
            }

            // Verificar se a célula a oeste pode ser contaminada
            if (coord.column > 0) {
                std::size_t west_key = key - 1;
                double tirage = pseudo_random(key*13427*13427*13427+m_time_step, m_time_step);
                double green_power = m_vegetation_map[west_key];
                double correction = power*log_factor(green_power);
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
                    local_updates.push_back({key, new_value});
                }
            } else {
                new_value = value >> 1; // Divide por 2
                local_updates.push_back({key, new_value});
            }
        }

        // Combinar resultados locais com os globais
        #pragma omp critical
        {
            new_fire_cells.insert(new_fire_cells.end(), local_new_fire.begin(), local_new_fire.end());
            updated_cells.insert(updated_cells.end(), local_updates.begin(), local_updates.end());
        }
    }

    // Aplicar as mudanças em uma única passagem
    // Primeiro as novas células em chamas
    for (const auto& key : new_fire_cells) {
        m_fire_map[key] = 255;
        m_fire_front[key] = 255;
    }

    // Depois as atualizações nas células existentes
    for (const auto& update : updated_cells) {
        std::size_t key = update.first;
        std::uint8_t new_value = update.second;
        
        m_fire_map[key] = new_value;
        if (new_value == 0) {
            m_fire_front.erase(key);
        } else {
            m_fire_front[key] = new_value;
        }
    }

    // Atualização da vegetação - com menos sincronização
    std::vector<std::size_t> vegetation_updates;
    vegetation_updates.reserve(fire_front_keys.size());

    #pragma omp parallel for
    for (size_t i = 0; i < fire_front_keys.size(); i++) {
        std::size_t key = fire_front_keys[i];
        if (m_vegetation_map[key] > 0) {
            #pragma omp atomic update
            m_vegetation_map[key] -= 1;
        }
    }

    m_time_step += 1;
    if (!m_fire_front.empty())
    {
        std::cout << "Fogo ainda ativo, continuando a simulação." << std::endl;
        return true;
    }
    else
    {
        std::cout << "Fogo extinto, encerrando a simulação." << std::endl;
        return false;
    }
    // return !m_fire_front.empty();
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