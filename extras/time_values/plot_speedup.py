import pandas as pd
import matplotlib.pyplot as plt

def read_time_data(filename):
    """
    Lê os dados de um arquivo `temps.txt` e retorna um DataFrame.
    """
    try:
        # Ler o arquivo, considerando que a primeira linha é o cabeçalho
        data = pd.read_csv(filename, delim_whitespace=True)
        return data
    except FileNotFoundError:
        print(f"Erreur: Fichier '{filename}' n'est pas trouvé")
        exit(1)
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier '{filename}': {e}")
        exit(1)

def adjust_timesteps(data):
    """
    Substitui os valores de TimeStep por uma sequência começando em 1,
    somente se o primeiro valor de TimeStep for diferente de 1.
    """
    if data["TimeStep"].iloc[0] != 1:  # Verifica se o primeiro valor de TimeStep é diferente de 1
        data["AdjustedTimeStep"] = range(1, len(data) + 1)  # Cria uma nova coluna com valores sequenciais
    else:
        data["AdjustedTimeStep"] = data["TimeStep"]  # Mantém os valores originais se o primeiro for 1
    return data

def filter_common_steps(sequential_data, parallel_data):
    """
    Filtra os dados para manter apenas os passos de tempo comuns entre os dois DataFrames.
    """
    # Encontrar a interseção dos passos de tempo
    common_steps = set(sequential_data["AdjustedTimeStep"]).intersection(set(parallel_data["AdjustedTimeStep"]))
    # Filtrar os dados para incluir apenas os passos de tempo comuns
    sequential_filtered = sequential_data[sequential_data["AdjustedTimeStep"].isin(common_steps)]
    parallel_filtered = parallel_data[parallel_data["AdjustedTimeStep"].isin(common_steps)]
    # Ordenar os dados pelo AdjustedTimeStep para garantir consistência
    sequential_filtered = sequential_filtered.sort_values(by="AdjustedTimeStep")
    parallel_filtered = parallel_filtered.sort_values(by="AdjustedTimeStep")
    return sequential_filtered, parallel_filtered

def calculate_speedup(sequential_times, parallel_times):
    """
    Calcula o speedup para cada passo de tempo.
    """
    speedups = []
    for seq, par in zip(sequential_times, parallel_times):
        if par == 0:
            print("Temps parallèle nul trouvé. Impossible de calculer la vitesse.")
            speedups.append(0)
        else:
            speedups.append(seq / par)
    return speedups

def plot_speedup(time_steps, speedups):
    """
    Gera um gráfico do speedup ao longo do tempo.
    """
    plt.figure(figsize=(10, 6))
    plt.plot(time_steps, speedups, label="Speedup", color="blue")  # Removido o marker="o"
    plt.title("Speedup au fil du temps")
    plt.xlabel("Time Step (adjusted)")
    plt.ylabel("Speedup")
    plt.grid(True)
    plt.legend()
    plt.savefig("speedup_plot.png")  # Salvar o gráfico como imagem
    plt.show()

def calculate_speedup(sequential_times, parallel_times):
    """
    Calcula o speedup para cada passo de tempo.
    """
    speedups = []
    for seq, par in zip(sequential_times, parallel_times):
        if par == 0:
            print("Tempo paralelo zero encontrado. Impossível calcular speedup.")
            speedups.append(0)
        else:
            speedups.append(seq / par)
    return speedups

def calculate_average_speedup(speedups):
    """
    Calcula o speedup médio.
    """
    if len(speedups) == 0:
        print("Nenhum speedup disponível para calcular a média.")
        return 0
    return sum(speedups) / len(speedups)

def main():
    # Nome dos arquivos de entrada
    sequential_file = "temps_seq_400.txt"
    parallel_file = "./400/temps_paral_4.txt"

    # Ler os dados dos arquivos
    sequential_data = read_time_data(sequential_file)
    parallel_data = read_time_data(parallel_file)

    # Ajustar os valores de TimeStep, se necessário
    sequential_data = adjust_timesteps(sequential_data)
    parallel_data = adjust_timesteps(parallel_data)

    # Filtrar apenas os passos de tempo comuns
    sequential_data, parallel_data = filter_common_steps(sequential_data, parallel_data)

    # Verificar se há dados suficientes
    if len(sequential_data) == 0 or len(parallel_data) == 0:
        print("Il n'y a pas de pas de temps commun entre les fichiers séquentiels et les fichiers parallèles.")
        exit(1)

    # Extrair os tempos de simulação
    sequential_times = sequential_data["Temps_Avancement"].values
    parallel_times = parallel_data["Temps_Avancement"].values
    time_steps = sequential_data["AdjustedTimeStep"].values  # Usar os passos ajustados

    # Calcular o speedup
    speedups = calculate_speedup(sequential_times, parallel_times)

    # Calcular o speedup médio
    average_speedup = calculate_average_speedup(speedups)

    # Exibir o resultado
    print(f"Speedup Médio: {average_speedup:.2f}")

    # Plotar o gráfico
    plot_speedup(time_steps, speedups)

if __name__ == "__main__":
    main()