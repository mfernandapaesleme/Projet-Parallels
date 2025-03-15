#include <stdexcept>
#include <cmath>
#include <iostream>
#include "model.hpp"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <vector>
#include <unordered_map>


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
        m_fire_map(t_discretization*t_discretization, 0u)
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
// --------------------------------------------------------------------------------------------------------------------
bool 
Model::update()
{
    
    std::vector<std::pair<std::size_t, unsigned char>> fire_front_vec(m_fire_front.begin(), m_fire_front.end());
    
    // Criar um container thread-local para acumular as atualizações do fire front
    tbb::combinable<std::unordered_map<std::size_t, unsigned char>> local_next_front(
        []() { return std::unordered_map<std::size_t, unsigned char>(); }
    );
    
    // Processa cada elemento do fire front em paralelo
    tbb::parallel_for(tbb::blocked_range<size_t>(0, fire_front_vec.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            // Obter o container local para esta thread
            auto &local_map = local_next_front.local();
            
            for (size_t i = r.begin(); i != r.end(); ++i) {
                auto f = fire_front_vec[i];
                LexicoIndices coord = get_lexicographic_from_index(f.first);
                double power = log_factor(f.second);
                
                // Processa a vizinhança: Sul
                if (coord.row < m_geometry - 1) {
                    double tirage = pseudo_random(f.first + m_time_step, m_time_step);
                    double green_power = m_vegetation_map[f.first + m_geometry];
                    double correction = power * log_factor(green_power);
                    if (tirage < alphaSouthNorth * p1 * correction) {
                        m_fire_map[f.first + m_geometry] = 255;
                        local_map[f.first + m_geometry] = 255;
                    }
                }
                
                // Norte
                if (coord.row > 0) {
                    double tirage = pseudo_random(f.first * 13427 + m_time_step, m_time_step);
                    double green_power = m_vegetation_map[f.first - m_geometry];
                    double correction = power * log_factor(green_power);
                    if (tirage < alphaNorthSouth * p1 * correction) {
                        m_fire_map[f.first - m_geometry] = 255;
                        local_map[f.first - m_geometry] = 255;
                    }
                }
                
                // Leste
                if (coord.column < m_geometry - 1) {
                    double tirage = pseudo_random(f.first * 13427 * 13427 + m_time_step, m_time_step);
                    double green_power = m_vegetation_map[f.first + 1];
                    double correction = power * log_factor(green_power);
                    if (tirage < alphaEastWest * p1 * correction) {
                        m_fire_map[f.first + 1] = 255;
                        local_map[f.first + 1] = 255;
                    }
                }
                
                // Oeste
                if (coord.column > 0) {
                    double tirage = pseudo_random(f.first * 13427 * 13427 * 13427 + m_time_step, m_time_step);
                    double green_power = m_vegetation_map[f.first - 1];
                    double correction = power * log_factor(green_power);
                    if (tirage < alphaWestEast * p1 * correction) {
                        m_fire_map[f.first - 1] = 255;
                        local_map[f.first - 1] = 255;
                    }
                }
                
                // Processa o enfraquecimento do fogo
                if (f.second == 255) {
                    double tirage = pseudo_random(f.first * 52513 + m_time_step, m_time_step);
                    if (tirage < p2) {
                        m_fire_map[f.first] >>= 1;
                        local_map[f.first] = f.second >> 1;
                    }
                } else {
                    m_fire_map[f.first] >>= 1;
                    local_map[f.first] = f.second >> 1;
                    // Se o fogo se extinguir (valor zero), não o adicionamos (ou removeremos depois da combinação)
                }
            }
        }
    );
    
    // Combina os resultados dos containers locais em um único container global para o fire front
    std::unordered_map<std::size_t, unsigned char> next_front;
    local_next_front.combine_each([&](const std::unordered_map<std::size_t, unsigned char>& local_map) {
        for (const auto &elem : local_map) {
            next_front[elem.first] = elem.second;
        }
    });
    
    // Atualiza o container global do fire front
    m_fire_front = std::move(next_front);
    
    // Atualiza a vegetação nas posições onde há fogo
    for(auto &f : m_fire_front) {
        if (m_vegetation_map[f.first] > 0)
            m_vegetation_map[f.first] -= 1;
    }
    
    m_time_step += 1;
    return !m_fire_front.empty();

}
// ====================================================================================================================
std::size_t   
Model::get_index_from_lexicographic_indices( LexicoIndices t_lexico_indices  ) const
{
    return t_lexico_indices.row*this->geometry() + t_lexico_indices.column;
}
// --------------------------------------------------------------------------------------------------------------------
auto 
Model::get_lexicographic_from_index( std::size_t t_global_index ) const -> LexicoIndices
{
    LexicoIndices ind_coords;
    ind_coords.row    = t_global_index/this->geometry();
    ind_coords.column = t_global_index%this->geometry();
    return ind_coords;
}
