import matplotlib.pyplot as plt
import numpy as np

# Função para calcular estatísticas sem outliers usando o método IQR
def calcular_sem_outliers(data):
    Q1 = np.percentile(data, 25)  
    Q3 = np.percentile(data, 75)  
    IQR = Q3 - Q1  
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return data[(data >= lower_bound) & (data <= upper_bound)]

# Ler os dados do arquivo
data = np.loadtxt("temps.txt", skiprows=1)

time_step = data[:, 0]
tempo_atualizacao = data[:, 1]  # Mantendo todos os dados para o gráfico
tempo_renderizacao = data[:, 2]
tempo_total = data[:, 3]

# Calcular estatísticas sem outliers
tempo_max_avanco_fogo = np.max(calcular_sem_outliers(tempo_atualizacao))
tempo_medio_display = np.mean(calcular_sem_outliers(tempo_renderizacao))

# Criar figuras separadas para cada métrica
plt.figure(figsize=(10, 7))

# Gráfico 1: Tempo de Avanço do Fogo
plt.subplot(2, 1, 1)
plt.plot(tempo_total, tempo_atualizacao, color="red")
plt.xlabel("Temps de la simulation")
plt.ylabel("Temp (s)")
plt.title("Temps pour l'avancement du feu")
plt.grid()

# Escrever o tempo máximo no gráfico
plt.text(
    0.05, 0.9, 
    f"Temps max (sans outliers): {tempo_max_avanco_fogo:.4f} s", 
    transform=plt.gca().transAxes, 
    fontsize=10, color="black"
)

# Gráfico 2: Tempo de Renderização
plt.subplot(2, 1, 2)
plt.plot(tempo_total, tempo_renderizacao, color="blue")
plt.xlabel("Temps de la simulation")
plt.ylabel("Temp (s)")
plt.title("Temps pour l'affichage")
plt.grid()

# Escrever o tempo médio no gráfico
plt.text(
    0.05, 0.9, 
    f"Temps moyen (sans outliers): {tempo_medio_display:.4f} s", 
    transform=plt.gca().transAxes, 
    fontsize=10, color="black"
)

# Ajustar espaçamento e salvar
plt.subplots_adjust(hspace=0.5)
plt.savefig("plot_temps.png")
plt.show()