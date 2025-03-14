#include <stdexcept>
#include <cmath>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <omp.h>
#include "model.hpp"

namespace
{
    double pseudo_random(std::size_t index, std::size_t time_step)
    {
        // Multiply index*((time_step+1) + 10000) removes
        std::uint_fast32_t xi = std::uint_fast32_t(index * (10000 + time_step + 1));
        std::uint_fast32_t r = (48271 * xi) % 2147483647;
        return r / 2147483646.;
    }

    double log_factor(std::uint8_t value)
    {
        return std::log(1.0 + value) / std::log(256.0);
    }
}

Model::Model(double t_length, unsigned t_discretization, std::array<double, 2> t_wind,
             LexicoIndices t_start_fire_position, double t_max_wind)
    : m_length(t_length),
      m_geometry(t_discretization),
      m_wind(t_wind),
      m_max_wind(t_max_wind),
      m_vegetation_map(t_discretization * t_discretization, 255u),
      m_fire_map(t_discretization * t_discretization, 0u),
      m_time_step(0)
{
    if (t_discretization == 0)
    {
        throw std::range_error("O número de células deve ser maior que zero.");
    }
    
    m_distance = m_length / double(m_geometry);
    auto index = get_index_from_lexicographic_indices(t_start_fire_position);
    m_fire_map[index] = 255u;
    m_fire_front[index] = 255u;

    double wind_speed = std::sqrt(t_wind[0] * t_wind[0] + t_wind[1] * t_wind[1]);
    p1 = (wind_speed < t_max_wind) ? (0.45279 + 0.000958 * wind_speed + 0.000036 * wind_speed * wind_speed)
                                   : (0.45279 + 0.000958 * t_max_wind + 0.000036 * t_max_wind * t_max_wind);
    p2 = 0.3;
}

bool Model::update()
{
    std::vector<std::size_t> fire_front_keys;
    fire_front_keys.reserve(m_fire_front.size());
    for (const auto &f : m_fire_front)
    {
        fire_front_keys.push_back(f.first);
    }

    if (fire_front_keys.empty())
        return false;

    std::vector<double> random_values(fire_front_keys.size());
    #pragma omp parallel for
    for (size_t i = 0; i < fire_front_keys.size(); ++i)
    {
        random_values[i] = pseudo_random(fire_front_keys[i], m_time_step);
    }

    std::unordered_map<std::size_t, std::uint8_t> next_front;
    std::vector<std::pair<std::size_t, std::uint8_t>> local_updates;
    
    #pragma omp parallel
    {
        std::vector<std::pair<std::size_t, std::uint8_t>> thread_updates;
        #pragma omp for
        for (size_t i = 0; i < fire_front_keys.size(); ++i)
        {
            std::size_t key = fire_front_keys[i];
            std::uint8_t value = m_fire_front[key];
            LexicoIndices coord = get_lexicographic_from_index(key);
            double power = log_factor(value);

            if (coord.row < m_geometry - 1)
            {
                double green_power = m_vegetation_map[key + m_geometry];
                if (random_values[i] < p1 * power * log_factor(green_power))
                {
                    thread_updates.emplace_back(key + m_geometry, 255u);
                }
            }
            if (coord.row > 0)
            {
                double green_power = m_vegetation_map[key - m_geometry];
                if (random_values[i] < p1 * power * log_factor(green_power))
                {
                    thread_updates.emplace_back(key - m_geometry, 255u);
                }
            }
            if (coord.column < m_geometry - 1)
            {
                double green_power = m_vegetation_map[key + 1];
                if (random_values[i] < p1 * power * log_factor(green_power))
                {
                    thread_updates.emplace_back(key + 1, 255u);
                }
            }
            if (coord.column > 0)
            {
                double green_power = m_vegetation_map[key - 1];
                if (random_values[i] < p1 * power * log_factor(green_power))
                {
                    thread_updates.emplace_back(key - 1, 255u);
                }
            }

            if (value == 255 && random_values[i] < p2)
                thread_updates.emplace_back(key, value >> 1);
            else
                thread_updates.emplace_back(key, value >> 1);
        }
        
        #pragma omp critical
        {
            local_updates.insert(local_updates.end(), thread_updates.begin(), thread_updates.end());
        }
    }

    for (const auto &[key, value] : local_updates)
    {
        m_fire_map[key] = value;
        if (value > 0)
        {
            next_front[key] = value;
        }
    }

    #pragma omp parallel for
    for (size_t i = 0; i < fire_front_keys.size(); ++i)
    {
        std::size_t key = fire_front_keys[i];
        if (m_vegetation_map[key] > 0)
        {
            #pragma omp atomic
            m_vegetation_map[key] -= 1;
        }
    }

    m_fire_front = next_front;
    m_time_step += 1;
    return !m_fire_front.empty();
}

std::size_t Model::get_index_from_lexicographic_indices(LexicoIndices t_lexico_indices) const
{
    return t_lexico_indices.row * this->geometry() + t_lexico_indices.column;
}

auto Model::get_lexicographic_from_index(std::size_t t_global_index) const -> LexicoIndices
{
    return {t_global_index / this->geometry(), t_global_index % this->geometry()};
}
