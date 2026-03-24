import torch
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import numpy as np
import os

class RadonDataset(Dataset):
    def __init__(self, csv_file, root_dir, transform=None):
        self.labels_frame = pd.read_csv(csv_file)
        self.root_dir = root_dir
        self.transform = transform
        
        # 1. Définition des colonnes de labels (Cible de l'IA)
        self.label_columns = ['vp_x', 'vp_y', 'vp_z']
        
        # 2. Identification dynamique du nombre de plans
        # On cherche toutes les colonnes qui finissent par '_x' (ex: n0_x, n1_x...)
        # en excluant 'vp_x'
        self.normal_cols = [c for c in self.labels_frame.columns if c.endswith('_x') and c != 'vp_x']
        self.nb_planes = len(self.normal_cols)
        
        # 3. Statistiques pour la normalisation Z-score des labels
        raw_labels = self.labels_frame[self.label_columns].values.astype('float32')
        self.labels_mean = np.mean(raw_labels, axis=0)
        self.labels_std = np.std(raw_labels, axis=0)
        self.labels_std[self.labels_std == 0] = 1e-8 

        print(f"📊 Dataset initialisé :")
        print(f"   Nombre de plans détectés : {self.nb_planes}")
        print(f"   Moyennes VP : {self.labels_mean}")
        print(f"   Écarts-types VP : {self.labels_std}")

    def __len__(self):
        return len(self.labels_frame)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()

        row = self.labels_frame.iloc[idx]

        # --- A. CHARGEMENT DES FEATURES (ENTRÉE) ---
        
        # 1. Lecture du binaire
        sample_id = int(row['id'])
        bin_path = os.path.join(self.root_dir, f"sample_{sample_id:06d}.bin")
        
        # Le binaire contient : nb_planes * (16*16*16) doubles
        # On utilise float64 car c'est le type 'double' du C++
        data = np.fromfile(bin_path, dtype=np.float64).astype(np.float32)
        
        # Reshape dynamique selon le nombre de plans détectés
        expected_size = self.nb_planes * (16**3)
        if data.size != expected_size:
            raise ValueError(f"Taille incorrecte pour {sample_id} : {data.size} au lieu de {expected_size}")
        
        data = data.reshape(self.nb_planes, 16, 16, 16)
        
        # Normalisation Min-Max locale par échantillon
        d_min, d_max = data.min(), data.max()
        if d_max > d_min:
            data = (data - d_min) / (d_max - d_min)
        
        # 2. Extraction des Normales de capture (Le "où")
        # On récupère les normales recentrées par Claire (nx, ny, nz pour chaque plan)
        normals = []
        for n in range(self.nb_planes):
            normals.extend([row[f'n{n}_x'], row[f'n{n}_y'], row[f'n{n}_z']])
        normals = np.array(normals, dtype='float32')

        # 3. Extraction de l'incrément (Le "zoom")
        increment = np.array([row['increment']], dtype='float32')

        # --- B. CHARGEMENT DES LABELS (CIBLE) ---
        raw_label = row[self.label_columns].values.astype('float32')
        normalized_label = (raw_label - self.labels_mean) / self.labels_std

        # --- C. PRÉPARATION DU DICTIONNAIRE ---
        # On renvoie un dictionnaire pour que le modèle sache qui est quoi
        return {
            'radon_cubes': torch.from_numpy(data),      # (NB_PLAN, 16, 16, 16)
            'normals': torch.from_numpy(normals),       # (NB_PLAN * 3)
            'increment': torch.from_numpy(increment),   # (1)
            'target_vp': torch.from_numpy(normalized_label) # (3)
        }

    def denormalize_labels(self, normalized_tensor):
        mean_tensor = torch.tensor(self.labels_mean, device=normalized_tensor.device)
        std_tensor = torch.tensor(self.labels_std, device=normalized_tensor.device)
        return (normalized_tensor * std_tensor) + mean_tensor

# --- TEST DU CHARGEMENT ---
if __name__ == "__main__":
    # Assure-toi que le chemin pointe vers ton nouveau dossier de data
    dataset = RadonDataset(csv_file='data/labels.csv', root_dir='data/')
    dataloader = DataLoader(dataset, batch_size=4, shuffle=True)

    try:
        batch = next(iter(dataloader))
        print("\n✅ Succès du chargement !")
        print(f"Dimensions Cubes : {batch['radon_cubes'].shape}") 
        print(f"Dimensions Normales : {batch['normals'].shape}")
        print(f"Valeur Incrément : {batch['increment'][0].item()}")
        print(f"Cible (VP normalisé) : {batch['target_vp'][0]}")
        
        # Test dénormalisation
        real_vp = dataset.denormalize_labels(batch['target_vp'][0])
        print(f"Vraie valeur VP : {real_vp}")
        
    except Exception as e:
        import traceback
        traceback.print_exc()