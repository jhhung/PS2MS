from base import BaseDataLoader
from torch.utils.data import Dataset
from scipy.sparse import load_npz
from utils.util import get_keep
import torch
import numpy as np

class SpectrumDataset(Dataset):
    def __init__(self, x_path, y_path=None, mw_path=None, keep=None, device='cpu', transform=None, max_mw=None):
        x = load_npz(x_path)
        self.x = torch.tensor(x.toarray(), dtype=torch.float32, device=device)
        self.mw = torch.tensor(np.load(mw_path)[..., None], dtype=torch.float32, device=device)
        if max_mw is None:
            self.max_mw = torch.max(self.mw)
        else:
            self.max_mw = max_mw
        self.mw /= self.max_mw
        self.device = device
        if y_path is not None:
            y = load_npz(y_path)
            if keep is None:
                keep = list(range(y.shape[1])) # keep all
            y = y[:, keep]
            pos_y = np.sum(y, axis=0) / y.shape[0]
            self.pos_weight = np.nan_to_num((np.ones_like(pos_y) - pos_y)/pos_y) # in case there is zero pos
            self.y = torch.tensor(y.toarray(), dtype=torch.float32, device=device)
        else:
            self.y = torch.zeros((x.shape[0],1))
        self.transform = transform
    
    def __len__(self):
        return self.x.shape[0]
    
    def __getitem__(self, idx):
        x = self.x[idx]
        if self.transform:
            x = self.transform(x)
        # return x, self.y[idx], self.mw[idx]
        return x, self.y[idx], self.mw[0] # ignore mw
    
    def get_pos_weight(self):
        return self.pos_weight if self.y is not None else None

class SpectrumDataloader(BaseDataLoader):
    def __init__(self, x_data_dir, y_data_dir, mw_dir, keep_dir, batch_size, shuffle, validation_split, num_workers=0, device='cpu', transform=None, max_mw=None):
        self.x_dir = x_data_dir
        self.y_dir = y_data_dir
        self.mw_dir = mw_dir
        self.keep_dir = keep_dir
        self.dataset = SpectrumDataset(self.x_dir, self.y_dir, self.mw_dir, get_keep(self.keep_dir), device=device, transform=transform, max_mw=max_mw)
        super().__init__(self.dataset, batch_size, shuffle, validation_split, num_workers)