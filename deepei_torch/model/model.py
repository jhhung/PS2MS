import torch
import torch.nn as nn
from base import BaseModel

class Simple(BaseModel):
    def __init__(self, input_dim, output_dim, layers, res_width=1024, activ=nn.ReLU(), dropout=0.5):
        super(Simple, self).__init__()
        self.flatten = nn.Flatten()
        self.activ = activ
        self.dropout = dropout
        modules = [
            nn.Linear(input_dim, res_width),
            nn.BatchNorm1d(res_width),
            self.activ,
            nn.Dropout(self.dropout)
        ]
        for _ in range(layers):
            modules.append(nn.Sequential(
                nn.Linear(res_width, res_width),
                nn.BatchNorm1d(res_width),
                self.activ,
                nn.Dropout(self.dropout)
            ))
        self.hidden_layers = nn.Sequential(*modules)
        self.output_fp = nn.Linear(res_width, output_dim)

    def forward(self, x):
        x = self.flatten(x)
        x = self.hidden_layers(x)
        fp_logits = self.output_fp(x)
        return fp_logits
