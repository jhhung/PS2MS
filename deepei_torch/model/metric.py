import torch


def accuracy(output, target):
    with torch.no_grad():
        pred = torch.sigmoid(output)
        assert pred.shape == target.shape
        correct = torch.sum((pred > 0.5) == (target > 0.5)).item()
    return correct / torch.numel(pred)
