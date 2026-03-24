import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.colors import LinearSegmentedColormap

def generate_custom_dark_fire_cube(filename, plan_idx=0, box_size=16):
    if not os.path.exists(filename):
        print(f"Erreur : {filename} introuvable.")
        return

    # 1. Chargement et Reshape
    data = np.fromfile(filename, dtype=np.float64)
    nb_plans = len(data) // (box_size**3)
    dataset = data.reshape((nb_plans, box_size, box_size, box_size))
    
    cube = dataset[plan_idx]
    cube_norm = (cube - cube.min()) / (cube.max() - cube.min() + 1e-8)

    # 2. Création de TOUS les points (16x16x16 = 4096 points)
    z, y, x = np.indices(cube.shape)
    x_f, y_f, z_f = x.flatten(), y.flatten(), z.flatten()
    v_f = cube_norm.flatten()

    # 3. CRÉATION DE LA COLORMAP SUR MESURE (Jaune -> Orange -> Rouge -> Gris -> Noir)
    # Définition des couleurs clés (0 = Noir, 1 = Jaune Brillant)
    # 'red_dark' et 'gray_dark' sont des couleurs hexadécimales personnalisées
    colors = [
        (0.0, "black"),          # Points vides (v=0)
        (0.15, "#2F2F2F"),      # Gris très foncé (bruit de fond très faible)
        (0.4, "#800000"),       # Rouge foncé (cœur du signal)
        (0.65, "red"),          # Rouge vif
        (0.85, "orange"),       # Orange
        (1.0, "yellow")         # Points brillants (v=1)
    ]
    # Création de la colormap segmentée linéaire
    custom_cmap = LinearSegmentedColormap.from_list("dark_fire", colors, N=256)

    # 4. Tracé 3D
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # --- OPACITÉ ET TAILLE ---
    # Pour un rendu "foncé", on augmente l'opacité minimale (0.4) 
    # et la taille minimale (15) pour bien voir tous les points.
    alphas = np.clip(v_f * 0.5 + 0.4, 0.4, 0.9)
    sizes = v_f * 150 + 15

    # Assignation des couleurs basées sur la colormap personnalisée et l'alpha
    rgba_colors = custom_cmap(v_f)
    rgba_colors[:, -1] = alphas # Remplacement canal Alpha

    scatter = ax.scatter(x_f, y_f, z_f, 
                        c=rgba_colors, 
                        s=sizes, 
                        edgecolors='none', 
                        antialiased=True)

    # 5. Esthétique Cube Parfait
    ax.set_box_aspect([1, 1, 1])
    b = box_size - 1
    ax.set_xlim(0, b); ax.set_ylim(0, b); ax.set_zlim(0, b)
    
    # Ticks de grille (0, 4, 8, 12, 15) pour montrer la structure de la boîte
    ax.set_xticks(range(0, box_size, 4))
    ax.set_yticks(range(0, box_size, 4))
    ax.set_zticks(range(0, box_size, 4))

    # Titre pour le poster
    ax.set_title(f"Dense Dark Radon Cube ($16^3 = 4096$ points)\nCustom Dark Fire Colormap - Plan {plan_idx}", 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Vue esthétique
    ax.view_init(elev=25, azim=45)
    
    # Nettoyage des parois pour un fond blanc pur
    ax.xaxis.pane.fill = ax.yaxis.pane.fill = ax.zaxis.pane.fill = False
    
    # Ajout d'une grille légère pour repère spatial
    ax.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    
    output_name = f"custom_dark_fire_cube_{plan_idx}.png"
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    print(f"✅ Grille dense avec colormap sur mesure sauvegardée : {output_name}")
    plt.show()

# --- EXÉCUTION ---
target_file = "data/sample_000000.bin" 
generate_custom_dark_fire_cube(target_file, plan_idx=0)