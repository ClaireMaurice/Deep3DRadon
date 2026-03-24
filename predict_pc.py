import torch
import numpy as np
import pandas as pd
from model_architecture import RadonCNN3D  # Importe ta classe de modèle

def inference_pc(bin_file_path, geometric_data, model_path, stats):
    """
    Effectue la prédiction du PC pour une nouvelle donnée.
    
    Args:
        bin_file_path (str): Chemin vers le fichier .bin (les 25 cubes Radon)
        geometric_data (dict): Dictionnaire contenant 'normals' (75) et 'increment' (1)
        model_path (str): Chemin vers le fichier radon_cnn_best.pth
        stats (dict): Moyennes et Écarts-types du dataset d'entraînement pour dénormaliser
    """
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    # 1. Chargement de l'architecture et des poids
    nb_planes = 25
    model = RadonCNN3D(nb_planes=nb_planes)
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.to(device)
    model.eval()  # Mode évaluation (désactive Dropout et fige BatchNorm)

    # 2. Préparation des cubes Radon (Entrée visuelle)
    # Lecture du binaire (nb_planes * 16*16*16 doubles)
    data = np.fromfile(bin_file_path, dtype=np.float64).astype(np.float32)
    data = data.reshape(1, nb_planes, 16, 16, 16) # Batch size de 1
    
    # Normalisation Min-Max locale (comme dans le dataset d'entraînement)
    d_min, d_max = data.min(), data.max()
    if d_max > d_min:
        data = (data - d_min) / (d_max - d_min)
    
    # 3. Préparation des paramètres géométriques (Entrée numérique)
    normals = torch.tensor(geometric_data['normals'], dtype=torch.float32).unsqueeze(0).to(device)
    increment = torch.tensor([geometric_data['increment']], dtype=torch.float32).unsqueeze(0).to(device)
    radon_cubes = torch.from_numpy(data).to(device)

    # 4. Inférence (Le "chemin inverse")
    with torch.no_grad(): # Pas de calcul de gradient pour gagner de la vitesse
        prediction_normalized = model(radon_cubes, normals, increment)
    
    # 5. Dénormalisation (Retour aux coordonnées réelles)
    # On repasse des valeurs Z-score vers les millimètres/pixels
    pred_np = prediction_normalized.cpu().numpy()[0]
    real_vp = (pred_np * stats['std']) + stats['mean']
    
    return real_vp

# --- EXEMPLE D'UTILISATION ---
if __name__ == "__main__":
    # Ces valeurs doivent provenir de ton analyse statistique sur le dossier /data initial
    # (Les valeurs imprimées par ton dataset.py lors de l'initialisation)
    training_stats = {
        'mean': np.array([0.5, 0.5, 0.6]), # Exemple de moyennes VP
        'std': np.array([0.003, 0.003, 0.005]) # Exemple d'écarts-types VP
    }

    # Exemple de données géométriques pour un nouvel échantillon
    # Remplacer par les vraies normales de ton fichier labels.csv ou ton code C++
    mock_geometric_data = {
        'normals': np.random.rand(75).astype('float32'), 
        'increment': 0.02
    }

    try:
        pc_corrigé = inference_pc(
            bin_file_path='data_v2/sample_000000.bin',
            geometric_data=mock_geometric_data,
            model_path='models/radon_cnn_best_v3.pth',
            stats=training_stats
        )
        
        print("\n✅ Inférence réussie !")
        print(f"Coordonnées PC prédites : X={pc_corrigé[0]:.6f}, Y={pc_corrigé[1]:.6f}, Z={pc_corrigé[2]:.6f}")
        
    except Exception as e:
        print(f"❌ Erreur lors de l'inférence : {e}")