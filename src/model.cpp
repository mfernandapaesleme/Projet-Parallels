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
        std::uint_fast32_t xi = std::uint_fast32_t(index * (time_step + 1));
        std::uint_fast32_t r = (48271 * xi) % 2147483647;
        return r / 2147483646.;
    }

    double log_factor(std::uint8_t value)
    {
        return std::log(1. + value) / std::log(256);
    }
}

Model::Model(double t_length, unsigned t_discretization, std::array<double, 2> t_wind,
             LexicoIndices t_start_fire_position, double t_max_wind)
    : m_length(t_length),
      m_distance(-1),
      m_geometry(t_discretization),
      m_wind(t_wind),
      m_wind_speed(std::sqrt(t_wind[0] * t_wind[0] + t_wind[1] * t_wind[1])),
      m_max_wind(t_max_wind),
      m_vegetation_map(t_discretization * t_discretization, 255u),
      m_fire_map(t_discretization * t_discretization, 0u),
      m_time_step(0)
{
    if (t_discretization == 0)
    {
        throw std::range_error("Le nombre de cases par direction doit être plus grand que zéro.");
    }
    m_distance = m_length / double(m_geometry);
    // Configuração do foco inicial
    auto index = get_index_from_lexicographic_indices(t_start_fire_position);
    std::cout << "Foco inicial configurado na posição: (" << t_start_fire_position.row << ", "
              << t_start_fire_position.column << ") com índice global: " << index << std::endl;

    m_fire_map[index] = 255u;
    m_fire_front[index] = 255u;

    constexpr double alpha0 = 4.52790762e-01;
    constexpr double alpha1 = 9.58264437e-04;
    constexpr double alpha2 = 3.61499382e-05;

    if (m_wind_speed < t_max_wind)
        p1 = alpha0 + alpha1 * m_wind_speed + alpha2 * (m_wind_speed * m_wind_speed);
    else
        p1 = alpha0 + alpha1 * t_max_wind + alpha2 * (t_max_wind * t_max_wind);
    p2 = 0.3;

    std::cout << "p1 (propagação): " << p1 << ", p2 (extinção): " << p2 << std::endl;

    if (m_wind[0] > 0)
    {
        alphaEastWest = std::abs(m_wind[0] / t_max_wind) + 1;
        alphaWestEast = 1. - std::abs(m_wind[0] / t_max_wind);
    }
    else
    {
        alphaWestEast = std::abs(m_wind[0] / t_max_wind) + 1;
        alphaEastWest = 1. - std::abs(m_wind[0] / t_max_wind);
    }

    if (m_wind[1] > 0)
    {
        alphaSouthNorth = std::abs(m_wind[1] / t_max_wind) + 1;
        alphaNorthSouth = 1. - std::abs(m_wind[1] / t_max_wind);
    }
    else
    {
        alphaNorthSouth = std::abs(m_wind[1] / t_max_wind) + 1;
        alphaSouthNorth = 1. - std::abs(m_wind[1] / t_max_wind);
    }
}

bool Model::update()
{
  std::cout << "Entrou na função update()" << std::endl;
  // Extrair todas as chaves do m_fire_front para um vetor
  std::vector<std::size_t> fire_front_keys;
  fire_front_keys.reserve(m_fire_front.size());
  for (const auto &f : m_fire_front)
  {
    fire_front_keys.push_back(f.first);
  }

  if (fire_front_keys.empty())
    {
        std::cout << "Foco inicial do incêndio não encontrado!" << std::endl;
        return false;
    }

  // Pré-calcular valores pseudo-aleatórios
  std::vector<double> random_values(fire_front_keys.size());
  #pragma omp parallel for
  for (size_t i = 0; i < fire_front_keys.size(); ++i)
  {
    random_values[i] = pseudo_random(fire_front_keys[i], m_time_step);
  }

  // Criar um novo front para a próxima iteração
  std::unordered_map<std::size_t, std::uint8_t> next_front;
  std::vector<std::pair<std::size_t, std::uint8_t>> local_updates;      

  // Processar cada célula em paralelo com OpenMP
  #pragma omp parallel
  {
    std::vector<std::pair<std::size_t, std::uint8_t>> thread_updates;
    #pragma omp for
    for (size_t i = 0; i < fire_front_keys.size(); ++i)
    {
      std::size_t key = fire_front_keys[i];
      std::uint8_t value = m_fire_front[key];

      // Recuperação da coordenada lexicográfica da case em fogo
      LexicoIndices coord = get_lexicographic_from_index(key);
      double power = log_factor(value);

      // Verificar células vizinhas
      if (coord.row < m_geometry - 1)
      {
        double green_power = m_vegetation_map[key + m_geometry];
        double correction = power * log_factor(green_power);
        if (random_values[i] < alphaSouthNorth * p1 * correction)
        {
          thread_updates.emplace_back(key + m_geometry, 255u);
        }
      }
      if (coord.row > 0)
      {
        double green_power = m_vegetation_map[key - m_geometry];
        double correction = power * log_factor(green_power);
        if (random_values[i] < alphaNorthSouth * p1 * correction)
        {
          thread_updates.emplace_back(key - m_geometry, 255u);
        }
      }
      if (coord.column < m_geometry - 1)
      {
        double green_power = m_vegetation_map[key + 1];
        double correction = power * log_factor(green_power);
        if (random_values[i] < alphaEastWest * p1 * correction)
        {
          thread_updates.emplace_back(key + 1, 255u);
        }
      }
      if (coord.column > 0)
      {
        double green_power = m_vegetation_map[key - 1];
        double correction = power * log_factor(green_power);
        if (random_values[i] < alphaWestEast * p1 * correction)
        {
          thread_updates.emplace_back(key - 1, 255u);
        }
      }
      // Atualizar estado da célula atual
      if (value == 255)
      {
        if (random_values[i] < p2)
        {
          thread_updates.emplace_back(key, value >> 1);
        }
      }
      else
      {
        thread_updates.emplace_back(key, value >> 1);
      }
    }

    // Fusão das atualizações locais
    #pragma omp critical
    {
      local_updates.insert(local_updates.end(), thread_updates.begin(), thread_updates.end());
    }
    std::cout << "Número de atualizações locais: " << local_updates.size() << std::endl;
  }

  // Aplicar as atualizações globais
  for (const auto &[key, value] : local_updates)
        {
            m_fire_map[key] = value;
            if (value > 0)
            {
                next_front[key] = value;
            }
        }
  std::cout << "Tamanho do próximo front: " << next_front.size() << std::endl;

  // Atualizar vegetação
  #pragma omp parallel for
  for (size_t i = 0; i < fire_front_keys.size(); ++i)
  {
    std::size_t key = fire_front_keys[i];
    if (m_vegetation_map[key] > 0)
    {
      #pragma omp atomic update
      m_vegetation_map[key] -= 1;
    }
    std::cout << "Densidade de vegetação após atualização: " << m_vegetation_map[key] << std::endl;
  }

  m_fire_front = next_front;
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
}

std::size_t Model::get_index_from_lexicographic_indices(LexicoIndices t_lexico_indices) const
{
    return t_lexico_indices.row * this->geometry() + t_lexico_indices.column;
}

auto Model::get_lexicographic_from_index(std::size_t t_global_index) const -> LexicoIndices
{
    LexicoIndices ind_coords;
    ind_coords.row = t_global_index / this->geometry();
    ind_coords.column = t_global_index % this->geometry();
    return ind_coords;
}