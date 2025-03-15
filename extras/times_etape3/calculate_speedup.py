import sys

def parse_log(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Pular a primeira linha (headers)
    data = lines[1:]
    
    update_times = []
    for line in data:
        parts = line.strip().split()
        if len(parts) >= 2:
            try:
                update_time = float(parts[1])
                update_times.append(update_time)
            except ValueError:
                continue
    if len(update_times) == 0:
        return None
    return sum(update_times) / len(update_times)

def main():
    # Tempo sequencial
    sequential_log = 'sequential.log'
    seq_update = parse_log(sequential_log)
    if seq_update is None:
        print("Erro: Não foi possível analisar o arquivo sequencial.")
        sys.exit(1)
    
    print(f'Tempo Update Sequencial Médio: {seq_update:.6f} segundos')
    
    # Tempos paralelos
    speedup_results = {}
    for threads in [2, 4, 8, 16]:
        parallel_log = f'parallel_{threads}.log'
        par_update = parse_log(parallel_log)
        if par_update is None:
            print(f'Erro: Não foi possível analisar o arquivo {parallel_log}.')
            continue
        speedup = seq_update / par_update
        speedup_results[threads] = speedup
        print(f'Tempo Update Paralelo Médio com {threads} threads: {par_update:.6f} segundos, Speedup: {speedup:.2f}')
    
    # Opcionalmente, salvar os resultados em um arquivo CSV
    with open('speedup_results.csv', 'w') as f:
        f.write('Threads,Speedup\n')
        for threads, sp in speedup_results.items():
            f.write(f'{threads},{sp:.2f}\n')

if __name__ == "__main__":
    main()