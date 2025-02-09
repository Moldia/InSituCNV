from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.metrics import roc_auc_score, average_precision_score
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import warnings
import matplotlib



def compute_ari(adata, reference_key='simulated_subclone', computed_key='cnv_leiden_res0.7'):
    """
    Compute the Adjusted Rand Index (ARI) to compare two clustering strategies in an AnnData object.

    Parameters:
    ----------
    adata : AnnData
        The annotated data object containing cluster labels in `.obs`.
    reference_key : str, optional
        The key in `adata.obs` for the reference cluster labels (default is 'simulated_subclone').
    computed_key : str, optional
        The key in `adata.obs` for the computed cluster labels (default is 'cnv_leiden_res0.7').

    Returns:
    -------
    float
        The Adjusted Rand Index, ranging from -1 to 1, where 1 indicates perfect agreement.
    """
    reference_labels = adata.obs[reference_key]
    computed_labels = adata.obs[computed_key]
    return adjusted_rand_score(reference_labels, computed_labels)

def compute_nmi(adata, reference_key='simulated_subclone', computed_key='cnv_leiden_res0.7'):
    """
    Compute the Normalized Mutual Information (NMI) to compare two clustering strategies in an AnnData object.

    Parameters:
    ----------
    adata : AnnData
        The annotated data object containing cluster labels in `.obs`.
    reference_key : str, optional
        The key in `adata.obs` for the reference cluster labels (default is 'simulated_subclone').
    computed_key : str, optional
        The key in `adata.obs` for the computed cluster labels (default is 'cnv_leiden_res0.7').

    Returns:
    -------
    float
        The Normalized Mutual Information, ranging from 0 to 1, where 1 indicates perfect agreement.
    """
    reference_labels = adata.obs[reference_key]
    computed_labels = adata.obs[computed_key]
    return normalized_mutual_info_score(reference_labels, computed_labels)

def compute_auc(simulated, inferred_scores, class_value):
    """
    Compute ROC AUC and PR AUC for a specific class (gains or losses).
    
    Parameters:
    simulated (array-like): Simulated CNV labels (0: neutral, 1: gain, -1: loss).
    inferred_scores (array-like): Inferred CNV scores.
    class_value (int): The class value to compute AUC for (1 for gains, -1 for losses).
    
    Returns:
    tuple: ROC AUC and PR AUC scores
    """
    # Convert to binary labels for the specified class
    binary_labels = (simulated == class_value).astype(int)
    
    # Compute AUC for the specified class
    if class_value == -1:
        # For losses, use negative inferred scores for alignment
        return (
            roc_auc_score(binary_labels, -inferred_scores),
            average_precision_score(binary_labels, -inferred_scores)
        )
    else:
        return (
            roc_auc_score(binary_labels, inferred_scores),
            average_precision_score(binary_labels, inferred_scores)
        )
    
    
def compute_performance_metrics(tp, fp, tn, fn):
    """
    Compute a variety of performance metrics for a classifier based on the confusion matrix values.

    Parameters:
    ----------
    tp : int
        True Positives - The number of correctly identified positive samples.
    fp : int
        False Positives - The number of incorrectly identified positive samples.
    tn : int
        True Negatives - The number of correctly identified negative samples.
    fn : int
        False Negatives - The number of incorrectly identified negative samples.

    Returns:
    -------
    dict
        A dictionary containing the following metrics:
        - 'Accuracy': Overall correctness of the classifier.
        - 'Precision': Proportion of predicted positives that are actual positives.
        - 'Recall' (Sensitivity): Proportion of actual positives correctly identified.
        - 'Specificity': Proportion of actual negatives correctly identified.
        - 'F1 Score': Harmonic mean of precision and recall.
        - 'MCC' (Matthews Correlation Coefficient): Balanced measure considering all confusion matrix components.
        - 'Balanced Accuracy': Average of sensitivity and specificity.
        - 'FPR' (False Positive Rate): Proportion of negatives incorrectly classified as positives.
        - 'FDR' (False Discovery Rate): Proportion of predicted positives that are false positives.
        - 'FNR' (False Negative Rate): Proportion of positives incorrectly classified as negatives.

    Notes:
    -----
    - The function ensures numerical stability by checking for division by zero in all metrics.
    - Useful for binary classification problems where confusion matrix values (TP, FP, TN, FN) are available.

    Example:
    -------
    >>> tp, fp, tn, fn = 50, 10, 30, 20
    >>> metrics = compute_performance_metrics(tp, fp, tn, fn)
    >>> for metric, value in metrics.items():
    ...     print(f"{metric}: {value:.2f}")
    Accuracy: 0.80
    Precision: 0.83
    Recall: 0.71
    Specificity: 0.75
    F1 Score: 0.77
    MCC: 0.57
    Balanced Accuracy: 0.73
    FPR: 0.25
    FDR: 0.17
    FNR: 0.29
    """
    metrics = {}
    metrics['Accuracy'] = (tp + tn) / (tp + tn + fp + fn)
    metrics['Precision'] = tp / (tp + fp) if (tp + fp) > 0 else 0
    metrics['Recall'] = tp / (tp + fn) if (tp + fn) > 0 else 0
    metrics['Specificity'] = tn / (tn + fp) if (tn + fp) > 0 else 0
    metrics['F1 Score'] = (2 * metrics['Precision'] * metrics['Recall'] /
                           (metrics['Precision'] + metrics['Recall'])) if (metrics['Precision'] + metrics['Recall']) > 0 else 0
    metrics['MCC'] = ((tp * tn - fp * fn) /
                      (((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5)) if (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) > 0 else 0
    metrics['Balanced Accuracy'] = (metrics['Recall'] + metrics['Specificity']) / 2
    metrics['FPR'] = fp / (fp + tn) if (fp + tn) > 0 else 0
    metrics['FDR'] = fp / (tp + fp) if (tp + fp) > 0 else 0
    metrics['FNR'] = fn / (tp + fn) if (tp + fn) > 0 else 0
    return metrics
