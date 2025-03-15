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
        print(f"Erro: Arquivo '{filename}' não encontrado.")
        exit(1)
    except Exception as e:
        print(f"Erro ao ler o arquivo '{filename}': {e}")
        exit(1)

def filter_common_steps(sequential_data, parallel_data):
    """
    Filtra os dados para manter apenas os passos de tempo comuns entre os dois DataFrames.
    """
    # Encontrar a interseção dos passos de tempo
    common_steps = set(sequential_data["TimeStep"]).intersection(set(parallel_data["TimeStep"]))
    # Filtrar os dados para incluir apenas os passos de tempo comuns
    sequential_filtered = sequential_data[sequential_data["TimeStep"].isin(common_steps)]
    parallel_filtered = parallel_data[parallel_data["TimeStep"].isin(common_steps)]
    # Ordenar os dados pelo TimeStep para garantir consistência
    sequential_filtered = sequential_filtered.sort_values(by="TimeStep")
    parallel_filtered = parallel_filtered.sort_values(by="TimeStep")
    return sequential_filtered, parallel_filtered

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

def plot_speedup(time_steps, speedups):
    """
    Gera um gráfico do speedup ao longo do tempo.
    """
    plt.figure(figsize=(10, 6))
    plt.plot(time_steps, speedups, label="Speedup", color="blue", marker="o")
    plt.title("Speedup ao Longo do Tempo")
    plt.xlabel("Passo de Tempo")
    plt.ylabel("Speedup")
    plt.grid(True)
    plt.legend()
    plt.savefig("speedup_plot.png")  # Salvar o gráfico como imagem
    plt.show()

def main():
    # Nome dos arquivos de entrada
    sequential_file = "temps_seq_400.txt"
    parallel_file = "./400/temps_paral_4_2.txt"

    # Ler os dados dos arquivos
    sequential_data = read_time_data(sequential_file)
    parallel_data = read_time_data(parallel_file)

    # Filtrar apenas os passos de tempo comuns
    sequential_data, parallel_data = filter_common_steps(sequential_data, parallel_data)

    # Verificar se há dados suficientes
    if len(sequential_data) == 0 or len(parallel_data) == 0:
        print("Não há passos de tempo comuns entre os arquivos sequencial e paralelo.")
        exit(1)

    # Extrair os tempos de simulação
    sequential_times = sequential_data["Temps_Avancement"].values
    parallel_times = parallel_data["Temps_Avancement"].values
    time_steps = sequential_data["Temps_Simulation"].values

    # Calcular o speedup
    speedups = calculate_speedup(sequential_times, parallel_times)

    # Plotar o gráfico
    plot_speedup(time_steps, speedups)

if __name__ == "__main__":
    main()