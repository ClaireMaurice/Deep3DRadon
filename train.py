import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, random_split
import os
from dataset import RadonDataset
from model_architecture import RadonCNN3D

# Hyperparamètres
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
EPOCHS = 100 # On peut monter un peu avec le scheduler
BATCH_SIZE = 16
LEARNING_RATE = 0.0005

# Chargement
full_dataset = RadonDataset(csv_file='data/labels.csv', root_dir='data/')
nb_planes = full_dataset.nb_planes # On récupère le nombre de plans détectés

train_size = int(0.8 * len(full_dataset))
val_size = len(full_dataset) - train_size
train_dataset, val_dataset = random_split(full_dataset, [train_size, val_size])

train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE, shuffle=False)

# Initialisation (On passe nb_planes au modèle)
model = RadonCNN3D(nb_planes=nb_planes).to(device)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)
scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=5)

print(f"🚀 Entraînement sur {nb_planes} plans | Device: {device}")

best_val_loss = float('inf')

for epoch in range(EPOCHS):
    model.train()
    running_loss = 0.0
    
    for batch in train_loader:
        # Extraction du dictionnaire
        cubes = batch['radon_cubes'].to(device)
        normals = batch['normals'].to(device)
        increment = batch['increment'].to(device)
        targets = batch['target_vp'].to(device)
        
        # Forward (Le modèle prend maintenant 3 arguments)
        outputs = model(cubes, normals, increment)
        loss = criterion(outputs, targets)
        
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
        running_loss += loss.item() * cubes.size(0)
    
    # Validation
    model.eval()
    val_loss = 0.0
    with torch.no_grad():
        for batch in val_loader:
            cubes, normals, increment, targets = batch['radon_cubes'].to(device), batch['normals'].to(device), batch['increment'].to(device), batch['target_vp'].to(device)
            outputs = model(cubes, normals, increment)
            val_loss += criterion(outputs, targets).item() * cubes.size(0)
            
    epoch_val_loss = val_loss / val_size
    scheduler.step(epoch_val_loss)

    if epoch_val_loss < best_val_loss:
        best_val_loss = epoch_val_loss
        torch.save(model.state_dict(), "models/radon_cnn_best.pth")
        status = "⭐ Saved!"
    else: status = ""

    print(f"Epoch {epoch+1:02d} | Loss: {running_loss/train_size:.6f} | Val: {epoch_val_loss:.6f} {status}")