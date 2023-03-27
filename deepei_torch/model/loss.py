from torch.nn import BCEWithLogitsLoss


def MultiLabelBCELoss(output, target, pos_weight=None):
    return BCEWithLogitsLoss(pos_weight=pos_weight)(output, target)