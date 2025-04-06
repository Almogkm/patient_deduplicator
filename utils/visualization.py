"""
Visualization utilities for the Patient Deduplicator package.
"""

import time
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

logger = logging.getLogger('patient_deduplication.visualization')

def plot_potential_duplicates(potential_duplicates, patient_profiles, output_dir=None, max_pairs=20):
    """
    Visualize the potential duplicate patient pairs with enhanced information.
    
    Parameters:
    -----------
    potential_duplicates : list
        List of dictionaries containing potential duplicate pairs
    patient_profiles : dict
        Dictionary of patient profiles keyed by uniquePatientKey
    output_dir : str or Path, optional
        Directory to save the visualization
    max_pairs : int
        Maximum number of pairs to visualize
        
    Returns:
    --------
    str
        Path to the saved visualization file
    """
    if not potential_duplicates:
        logger.warning("No potential duplicates to visualize")
        return None
        
    # Set output directory
    if output_dir:
        output_dir = Path(output_dir)
    else:
        output_dir = Path.cwd() / "output"
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Prepare data for visualization
    pairs_to_viz = min(max_pairs, len(potential_duplicates))
    top_pairs = potential_duplicates[:pairs_to_viz]
    
    # Create labels with more detailed information
    pair_labels = []
    similarity_scores = []
    match_types = []
    
    for p in top_pairs:
        profile1 = patient_profiles[p['patient1']]
        profile2 = patient_profiles[p['patient2']]
        
        # Get patient IDs
        patient1_id = profile1.get('patient_id', p['patient1'])
        patient2_id = profile2.get('patient_id', p['patient2'])
        
        # Get sex/gender if available
        sex1 = profile1.get('SEX', profile1.get('sex', profile1.get('GENDER', profile1.get('gender', ''))))
        sex2 = profile2.get('SEX', profile2.get('sex', profile2.get('GENDER', profile2.get('gender', ''))))
        
        # Get age if available
        age1 = profile1.get('AGE', profile1.get('age', ''))
        age2 = profile2.get('AGE', profile2.get('age', ''))
        
        # Create label
        label = f"{patient1_id[:10]}:{p['study1']}"
        if sex1 or age1:
            demo1 = []
            if sex1:
                demo1.append(str(sex1))
            if age1:
                demo1.append(f"Age:{age1}")
            label += f" ({', '.join(demo1)})"
        
        label += f" - {patient2_id[:10]}:{p['study2']}"
        if sex2 or age2:
            demo2 = []
            if sex2:
                demo2.append(str(sex2))
            if age2:
                demo2.append(f"Age:{age2}")
            label += f" ({', '.join(demo2)})"
        
        pair_labels.append(label)
        similarity_scores.append(p['similarity'])
        
        # Determine match type
        details = p.get('details', {})
        if isinstance(details, dict):
            if details.get('clinical') and details.get('mutation') and 'error' not in details.get('mutation', {}):
                match_types.append('Combined')
            elif details.get('mutation') and 'error' not in details.get('mutation', {}):
                match_types.append('Mutation only')
            elif details.get('clinical'):
                match_types.append('Clinical only')
            else:
                match_types.append('Unknown')
        else:
            match_types.append('Unknown')
    
    # Create visualization - horizontal bar chart with color coding
    plt.figure(figsize=(14, 10))
    
    # Define colors for different match types
    colors = {'Combined': 'darkblue', 'Mutation only': 'darkorange', 
              'Clinical only': 'darkgreen', 'Unknown': 'gray'}
    
    # Create bars with colors based on match type
    bars = plt.barh(pair_labels, similarity_scores, color=[colors[t] for t in match_types])
    
    # Add match type labels to the bars
    for i, (score, match_type) in enumerate(zip(similarity_scores, match_types)):
        plt.text(score + 0.02, i, match_type, va='center', fontsize=8)
    
    plt.xlabel('Similarity Score')
    plt.ylabel('Patient Pairs')
    plt.title('Top Potential Duplicate Patient Pairs')
    plt.xlim(0, 1.15)  # Extended to fit labels
    
    # Add a legend
    legend_handles = [plt.Rectangle((0, 0), 1, 1, color=colors[t]) for t in set(match_types)]
    plt.legend(legend_handles, set(match_types), loc='lower right')
    
    plt.tight_layout()
    
    # Save figure
    timestamp = time.strftime('%Y%m%d-%H%M%S')
    viz_path = output_dir / f'duplicate_pairs_visualization_{timestamp}.png'
    plt.savefig(viz_path)
    logger.info(f"Saved visualization to {viz_path}")
    
    plt.close()
    return str(viz_path)

def plot_similarity_distribution(potential_duplicates, output_dir=None):
    """
    Plot the distribution of similarity scores for potential duplicates.
    
    Parameters:
    -----------
    potential_duplicates : list
        List of dictionaries containing potential duplicate pairs
    output_dir : str or Path, optional
        Directory to save the visualization
        
    Returns:
    --------
    str
        Path to the saved visualization file
    """
    if not potential_duplicates:
        logger.warning("No potential duplicates to visualize")
        return None
        
    # Set output directory
    if output_dir:
        output_dir = Path(output_dir)
    else:
        output_dir = Path.cwd() / "output"
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Extract similarity scores
    scores = [p['similarity'] for p in potential_duplicates]
    
    # Create histogram
    plt.figure(figsize=(10, 6))
    plt.hist(scores, bins=20, color='skyblue', edgecolor='black')
    plt.axvline(x=0.7, color='red', linestyle='--', label='Threshold (0.7)')
    
    plt.xlabel('Similarity Score')
    plt.ylabel('Number of Pairs')
    plt.title('Distribution of Similarity Scores for Potential Duplicates')
    plt.legend()
    
    plt.grid(axis='y', alpha=0.75)
    plt.tight_layout()
    
    # Save figure
    timestamp = time.strftime('%Y%m%d-%H%M%S')
    hist_path = output_dir / f'similarity_distribution_{timestamp}.png'
    plt.savefig(hist_path)
    logger.info(f"Saved similarity distribution to {hist_path}")
    
    plt.close()
    return str(hist_path)

def plot_study_heatmap(potential_duplicates, output_dir=None):
    """
    Create a heatmap showing the number of potential duplicates between studies.
    
    Parameters:
    -----------
    potential_duplicates : list
        List of dictionaries containing potential duplicate pairs
    output_dir : str or Path, optional
        Directory to save the visualization
        
    Returns:
    --------
    str
        Path to the saved visualization file
    """
    if not potential_duplicates:
        logger.warning("No potential duplicates to visualize")
        return None
        
    # Set output directory
    if output_dir:
        output_dir = Path(output_dir)
    else:
        output_dir = Path.cwd() / "output"
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Extract study pairs
    study_pairs = [(p['study1'], p['study2']) for p in potential_duplicates]
    
    # Count occurrences of each study pair
    study_counts = {}
    for s1, s2 in study_pairs:
        # Ensure consistent ordering (s1 <= s2)
        if s1 > s2:
            s1, s2 = s2, s1
        
        pair = (s1, s2)
        study_counts[pair] = study_counts.get(pair, 0) + 1
    
    # Get unique studies
    all_studies = set()
    for s1, s2 in study_counts.keys():
        all_studies.add(s1)
        all_studies.add(s2)
    
    all_studies = sorted(all_studies)
    n_studies = len(all_studies)
    
    # Create matrix for heatmap
    heatmap_matrix = np.zeros((n_studies, n_studies))
    
    # Study to index mapping
    study_to_idx = {study: i for i, study in enumerate(all_studies)}
    
    # Fill matrix
    for (s1, s2), count in study_counts.items():
        i, j = study_to_idx[s1], study_to_idx[s2]
        heatmap_matrix[i, j] = count
        heatmap_matrix[j, i] = count  # Make it symmetric
    
    # Create heatmap
    plt.figure(figsize=(12, 10))
    plt.imshow(heatmap_matrix, cmap='YlOrRd')
    
    # Add text annotations
    for i in range(n_studies):
        for j in range(n_studies):
            if heatmap_matrix[i, j] > 0:
                plt.text(j, i, int(heatmap_matrix[i, j]), ha='center', va='center', 
                         color='black' if heatmap_matrix[i, j] < 10 else 'white')
    
    # Add labels
    plt.xticks(range(n_studies), all_studies, rotation=90)
    plt.yticks(range(n_studies), all_studies)
    
    plt.title('Potential Duplicates Between Studies')
    plt.colorbar(label='Number of Potential Duplicates')
    
    plt.tight_layout()
    
    # Save figure
    timestamp = time.strftime('%Y%m%d-%H%M%S')
    heatmap_path = output_dir / f'study_heatmap_{timestamp}.png'
    plt.savefig(heatmap_path)
    logger.info(f"Saved study heatmap to {heatmap_path}")
    
    plt.close()
    return str(heatmap_path)