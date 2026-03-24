import torch
import torch.nn as nn
import torch.nn.functional as F

class RadonCNN3D(nn.Module):
    def __init__(self, nb_planes=25):
        super(RadonCNN3D, self).__init__()
        
        # --- BRANCHE 1 : CNN 3D (Analyse des images Radon) ---
        self.conv1 = nn.Conv3d(in_channels=nb_planes, out_channels=64, kernel_size=3, padding=1)
        self.bn1 = nn.BatchNorm3d(64)
        self.pool1 = nn.MaxPool3d(kernel_size=2) 
        
        self.conv2 = nn.Conv3d(64, 128, kernel_size=3, padding=1)
        self.bn2 = nn.BatchNorm3d(128)
        self.pool2 = nn.MaxPool3d(kernel_size=2) 
        
        self.conv3 = nn.Conv3d(128, 256, kernel_size=3, padding=1)
        self.bn3 = nn.BatchNorm3d(256)
        self.pool3 = nn.MaxPool3d(kernel_size=2) 

        # --- BRANCHE 2 : MLP CONTEXTUEL (Analyse Normales + Incrément) ---
        # Entrée : (nb_planes * 3) normales + 1 incrément
        context_input_dim = (nb_planes * 3) + 1
        self.context_mlp = nn.Sequential(
            nn.Linear(context_input_dim, 64),
            nn.LeakyReLU(0.1),
            nn.Linear(64, 64),
            nn.LeakyReLU(0.1)
        )

        # --- FUSION ET DÉCISION ---
        # Sortie CNN (256*2*2*2 = 2048) + Sortie MLP (64)
        self.fc1 = nn.Linear(2048 + 64, 512) 
        self.dropout1 = nn.Dropout(0.4)
        
        self.fc2 = nn.Linear(512, 128)
        self.dropout2 = nn.Dropout(0.3)
        
        self.fc_out = nn.Linear(128, 3) 

    def forward(self, cubes, normals, increment):
        # 1. Traitement des cubes
        x1 = self.pool1(F.leaky_relu(self.bn1(self.conv1(cubes)), 0.1))
        x1 = self.pool2(F.leaky_relu(self.bn2(self.conv2(x1)), 0.1))
        x1 = self.pool3(F.leaky_relu(self.bn3(self.conv3(x1)), 0.1))
        x1 = x1.view(x1.size(0), -1) # Flatten (batch, 2048)
        
        # 2. Traitement du contexte (Normales + Incrément)
        # On concatène les normales et l'incrément
        context = torch.cat((normals, increment), dim=1) 
        x2 = self.context_mlp(context) # (batch, 64)
        
        # 3. Fusion
        combined = torch.cat((x1, x2), dim=1) # (batch, 2112)
        
        # 4. Décision finale
        x = F.leaky_relu(self.fc1(combined), 0.1)
        x = self.dropout1(x)
        x = F.leaky_relu(self.fc2(x), 0.1)
        x = self.dropout2(x)
        
        return self.fc_out(x)