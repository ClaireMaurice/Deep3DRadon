import torch
import pandas as pd
import numpy as np
import os
from torch.utils.data import DataLoader
from dataset import RadonDataset 
from model_architecture import RadonCNN3D 

def evaluate():
    # 1. Configuration
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    csv_file = 'data/labels.csv'
    root_dir = 'data/'
    # On charge le "Best" modèle (celui qui a la meilleure Val Loss)
    model_path = 'models/radon_cnn_best.pth' 
    
    # 2. Chargement du Dataset
    dataset = RadonDataset(csv_file=csv_file, root_dir=root_dir)
    test_loader = DataLoader(dataset, batch_size=5, shuffle=True)

    # 3. Initialisation et Chargement du Modèle
    model = RadonCNN3D().to(device)
    if not os.path.exists(model_path):
        print(f"❌ Erreur : Le fichier {model_path} est introuvable.")
        print("Assure-toi d'avoir lancé train.py au moins une fois.")
        return

    model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval() 

    print(f"✅ Modèle chargé sur {device}. Analyse de 5 échantillons aléatoires (Dénormalisés)...\n")
    print(f"{'ID':<12} | {'Réel (x, y, z)':<25} | {'Prédit (x, y, z)':<25} | {'Erreur L2':<10}")
    print("-" * 90)

    with torch.no_grad(): 
        # On récupère un batch
        inputs, targets_norm = next(iter(test_loader))
        inputs, targets_norm = inputs.to(device), targets_norm.to(device)
        
        # Prédiction (en sortie normalisée)
        outputs_norm = model(inputs)

        # --- ÉTAPE CRUCIALE : Dénormalisation ---
        # On repasse des unités statistiques aux unités réelles (ex: mm)
        targets_real = dataset.denormalize_labels(targets_norm)
        outputs_real = dataset.denormalize_labels(outputs_norm)

        # Calcul de l'erreur sur les valeurs réelles
        for i in range(len(targets_real)):
            real = targets_real[i].cpu().numpy()
            pred = outputs_real[i].cpu().numpy()
            
            # Distance Euclidienne en unités réelles
            error = np.linalg.norm(real - pred) 
            
            # Formatage des chaînes de caractères
            real_str = f"[{real[0]:.4f}, {real[1]:.4f}, {real[2]:.4f}]"
            pred_str = f"[{pred[0]:.4f}, {pred[1]:.4f}, {pred[2]:.4f}]"
            
            print(f"Sample {i:<5} | {real_str:<25} | {pred_str:<25} | {error:.6f}")

    print("-" * 90)
    print("\nNote : L'erreur L2 est maintenant exprimée dans l'unité de tes coordonnées CSV.")

if __name__ == "__main__":
    evaluate()