import seaborn as sns
import pandas as pd
import numpy as np
from scipy.sparse import load_npz
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score

def bits2array(bits):
    return np.array(list(map(int, bits)), dtype=np.int64)

def plot_Score(ans, pred_fp_deepei, pred_fp_deeptorch, filename, target):
    fp_len = pred_fp_deepei.shape[1]
    metrics = [(accuracy_score, 'Accuracy'), (precision_score, 'Precision'), (recall_score, 'Recall'), (f1_score, 'F1-Score')]
    impl = [('DeepEI', pred_fp_deepei), ('DeepEI torch', pred_fp_deeptorch)]

    zd = {'zero_division':0}
    placeholder = {'sample_weight':None}

    scores = []
    metric_names = []
    impl_names = []
    
    for impl_name, impl_pred in impl:
        for metr, metr_name in metrics:
            scores += [metr(ans[:, i], impl_pred[:, i], **zd if metr_name != metrics[0][1] else placeholder) for i in range(fp_len)]
            metric_names += [metr_name] * fp_len
        impl_names += [impl_name] * (fp_len * len(metrics))
        
    plot_df = pd.DataFrame({
        'Score': scores,
        'Metric': metric_names,
        'Impl': impl_names,
    })
    fig, ax = plt.subplots()
    ax.grid(axis='y')
    ax.set_axisbelow(True)
    ax.set_title('Metric score of DeepEI and Our impl on ' + target)
    sns.violinplot(data=plot_df, x='Metric', y='Score', hue="Impl", split=True, cut=0., inner='box', ax=ax)
    ax.set_ylim(0., 1.)
    ax.set_yticks(np.linspace(0, 1, 11))
    fig.tight_layout()
    fig.savefig(filename)

def plot_Score_diff(ans, pred_fp_deepei, pred_fp_deeptorch, filename, target):
    fp_len = pred_fp_deepei.shape[1]
    metrics = [(accuracy_score, 'Accuracy'), (precision_score, 'Precision'), (recall_score, 'Recall'), (f1_score, 'F1-Score')]
    impl = {'DeepEI': pred_fp_deepei, 'DeepEI torch': pred_fp_deeptorch}

    zd = {'zero_division':0}
    placeholder = {'sample_weight':None}

    scores = []
    metric_names = []
    
    for metr, metr_name in metrics:
        scores += [
            metr(ans[:, i], impl['DeepEI torch'][:, i],  **zd if metr_name != metrics[0][1] else placeholder) -
            metr(ans[:, i], impl['DeepEI'][:, i],        **zd if metr_name != metrics[0][1] else placeholder)
                for i in range(fp_len)]
        metric_names += [metr_name] * fp_len
        
    plot_df = pd.DataFrame({
        'Diff (DeepEI torch - DeepEI)': scores,
        'Metric': metric_names,
    })
    fig, ax = plt.subplots()
    ax.grid(axis='y')
    ax.set_axisbelow(True)
    ax.set_title('Metric score difference between DeepEI and Our impl on ' + target)
    sns.violinplot(data=plot_df, x='Metric', y='Diff (DeepEI torch - DeepEI)', cut=0., ax=ax)
    fig.tight_layout()
    fig.savefig(filename)

def plot_Confidence(ans, pred_fp_deepei, pred_fp_deeptorch, filename, target):

    pass

def main(target):
    pred_df_deepei = pd.read_csv(f'/mnt/ec/mammoth/blender/spectrum/deepei/pred/nist_only/{target}_pred.csv')
    pred_df_deeptorch = pd.read_csv(f'/mnt/ec/mammoth/blender/spectrum/deepei_torch/pred/nist_only/{target}_pred.csv')
    ans = load_npz(f'/mnt/ec/ness/blender/repos/DeepEI/DeepEI/data/paper/{target}_fingerprint.npz').toarray()
    pred_fp_deepei = np.stack(list(map(bits2array, pred_df_deepei['fingerprint'])))
    pred_fp_deeptorch = np.stack(list(map(bits2array, pred_df_deeptorch['fingerprint'])))
    print("data loaded")
    sample_count = pred_fp_deepei.shape[0]
    

    plot_Score(ans, pred_fp_deepei, pred_fp_deeptorch, f'img/{target}_Score_distr_pos_wise.png', target)

    plot_Score_diff(ans, pred_fp_deepei, pred_fp_deeptorch, f'img/{target}_Score_diff_pos_wise.png', target)
    return

if __name__ == "__main__":
    main('nist_test')
    main('swg_test')
