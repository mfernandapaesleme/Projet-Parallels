# Analyse des résultats
echo "=== Résultats d'accélération ==="
echo "Threads | Temps total | Accélération | Temps calcul | Accélération calcul"
echo "--------------------------------------------------------------------"

# Extraire les temps de l'exécution séquentielle
SEQ_TOTAL=$(grep "Temps total moyen" sequential.log | awk '{print $6}')
SEQ_COMPUTE=$(grep "Temps de calcul moyen" sequential.log | awk '{print $7}')

echo "1 | $SEQ_TOTAL | 1.00 | $SEQ_COMPUTE | 1.00"

# Calculer les accélérations pour les exécutions parallèles
for threads in 2 4 8 16
do
    if [ -f parallel_${threads}.log ]; then
        TOTAL=$(grep "Temps total moyen" parallel_${threads}.log | awk '{print $6}')
        COMPUTE=$(grep "Temps de calcul moyen" parallel_${threads}.log | awk '{print $7}')
        
        # Calculer les accélérations
        if [ ! -z "$TOTAL" ] && [ ! -z "$COMPUTE" ]; then
            SPEEDUP_TOTAL=$(echo "scale=2; $SEQ_TOTAL / $TOTAL" | bc)
            SPEEDUP_COMPUTE=$(echo "scale=2; $SEQ_COMPUTE / $COMPUTE" | bc)
            echo "$threads | $TOTAL | $SPEEDUP_TOTAL | $COMPUTE | $SPEEDUP_COMPUTE"
        else
            echo "$threads | Erreur dans les mesures"
        fi
    else
        echo "$threads | Fichier de log manquant"
    fi
done