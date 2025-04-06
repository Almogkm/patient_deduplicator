"""
Similarity computation utilities for the Patient Deduplicator package.
"""

import logging
import pandas as pd

logger = logging.getLogger('patient_deduplication.similarity')

def compute_clinical_similarity(profile1, profile2, min_matches=2):
    """
    Compute similarity between two patients based on their clinical profiles.
    
    Parameters:
    -----------
    profile1 : dict
        First patient's profile
    profile2 : dict
        Second patient's profile
    min_matches : int
        Minimum number of matching fields required for a valid comparison
        
    Returns:
    --------
    float
        Clinical similarity score between 0 and 1
    dict
        Detailed breakdown of similarity components
    """
    # Helper function to replace if value is NA-like with None
    def replace_unknown_value(val):
        unknown_variants = {'unknown', 'unkown', 'u', 'n/a', 'na', 'nan', ''}
        if pd.isna(val) or str(val).strip().lower() in unknown_variants:
            return None
        return val

    # Compare demographic fields if they exist in both profiles
    demographic_fields = ['AGE', 'SEX', 'GENDER', 'RACE', 'ETHNICITY']
    
    # Track non-empty fields and their values
    non_empty_fields = []
    field_values = {}
    
    # First pass: collect all valid demographic fields
    for field in demographic_fields:
        field_lower = field.lower()
        
        # Skip if field doesn't exist in either profile
        if field not in profile1 and field_lower not in profile1:
            continue
        if field not in profile2 and field_lower not in profile2:
            continue
        
        # Get field values, handling case variations
        value1 = profile1.get(field, profile1.get(field_lower))
        value2 = profile2.get(field, profile2.get(field_lower))
        
        # Skip if either value is NA
        value1 = replace_unknown_value(value1)
        value2 = replace_unknown_value(value2)
        if value1 is None or value2 is None:
            continue
        
        non_empty_fields.append(field)
        field_values[field] = (value1, value2)
    
    # Early exit if we don't have enough valid fields for clinical comparison
    if len(non_empty_fields) < min_matches:
        logger.debug(f"Insufficient non-empty fields ({len(non_empty_fields)}). Minimum required: {min_matches}")
        return 0.0, {'error': 'Insufficient clinical data'}
    
    # Compare exactly min_matches number of fields
    fields_to_compare = non_empty_fields[:min_matches]
    clinical_matches = 0
    clinical_details = {}
    
    for field in fields_to_compare:
        value1, value2 = field_values[field]
        
        # Exact match comparison for all fields
        if str(value1).upper() == str(value2).upper():
            clinical_matches += 1
            clinical_details[field] = 1.0
        else:
            clinical_details[field] = 0.0
    
    clinical_score = clinical_matches / min_matches
    
    details = {
        'matches': clinical_matches,
        'total_compared': min_matches,
        'non_empty_fields': len(non_empty_fields),
        'field_details': clinical_details
    }
    
    return clinical_score, details

def compute_mutation_similarity(mut_profile1, mut_profile2):
    """
    Compute similarity between two patients based on their mutation profiles.
    
    Parameters:
    -----------
    mut_profile1 : dict
        First patient's mutation profile
    mut_profile2 : dict
        Second patient's mutation profile
        
    Returns:
    --------
    float
        Mutation similarity score between 0 and 1
    dict
        Detailed breakdown of similarity components
    """
    # Get mutation signatures and genes
    sig1 = mut_profile1.get('mutation_signatures', set())
    sig2 = mut_profile2.get('mutation_signatures', set())
    genes1 = mut_profile1.get('mutated_genes', set())
    genes2 = mut_profile2.get('mutated_genes', set())
    
    # Error checks
    if not sig1 or not sig2:
        return 0.0, {'error': 'Empty mutation signatures'}
    
    # Compute Jaccard similarity for mutation signatures
    jaccard_sig = len(sig1.intersection(sig2)) / len(sig1.union(sig2)) if sig1 and sig2 else 0
    
    # Compute Jaccard similarity for mutated genes
    jaccard_genes = len(genes1.intersection(genes2)) / len(genes1.union(genes2)) if genes1 and genes2 else 0
    
    # Overall mutation score - weighted combination
    mutation_score = 0.7 * jaccard_sig + 0.3 * jaccard_genes
    
    # Mutation details
    details = {
        'mutation_signature_similarity': jaccard_sig,
        'mutated_genes_similarity': jaccard_genes,
        'common_mutations': len(sig1.intersection(sig2)),
        'total_unique_mutations': len(sig1.union(sig2)),
        'common_genes': len(genes1.intersection(genes2)),
        'total_unique_genes': len(genes1.union(genes2))
    }
    
    return mutation_score, details

def compute_combined_similarity(clinical_score, clinical_details, mutation_score, mutation_details):
    """
    Compute a combined similarity score using clinical and mutation data.
    
    Parameters:
    -----------
    clinical_score : float
        Clinical similarity score
    clinical_details : dict
        Details of clinical similarity
    mutation_score : float
        Mutation similarity score
    mutation_details : dict
        Details of mutation similarity
        
    Returns:
    --------
    float
        Combined similarity score between 0 and 1
    dict
        Detailed breakdown of similarity components
    """
    # Weight based on data availability
    mutation_weight = 0.6 if mutation_score is not None else 0
    clinical_weight = 0.4 if clinical_score is not None else 0
    
    # Adjust weights if only one type of data is available
    if mutation_weight == 0 and clinical_weight > 0:
        clinical_weight = 1.0
    elif clinical_weight == 0 and mutation_weight > 0:
        mutation_weight = 1.0
    
    # Calculate combined score
    if mutation_weight + clinical_weight > 0:
        combined_score = (
            (mutation_score * mutation_weight + 
             clinical_score * clinical_weight) / 
            (mutation_weight + clinical_weight)
        )
    else:
        combined_score = 0.0
    
    # Combined details
    details = {
        'clinical': clinical_details,
        'mutation': mutation_details,
        'weights': {
            'clinical': clinical_weight,
            'mutation': mutation_weight
        }
    }
    
    return combined_score, details

def find_potential_duplicates(patient_profiles, mutation_profiles=None, 
                             similarity_threshold=0.7, clinical_threshold=0.8,
                             min_matches=2, method='combined', max_comparisons=None):
    """
    Find potential duplicate patients across different studies.
    
    Parameters:
    -----------
    patient_profiles : dict
        Dictionary of patient profiles keyed by uniquePatientKey
    mutation_profiles : dict, optional
        Dictionary of mutation profiles keyed by uniquePatientKey
    similarity_threshold : float
        Threshold for considering patients as potential duplicates (0-1)
    clinical_threshold : float
        Higher threshold for clinical-only matches (0-1)
    min_matches : int
        Minimum number of matching fields required for a valid comparison
    method : str
        Similarity method to use ('combined', 'mutations_only', 'clinical_only')
    max_comparisons : int, optional
        Maximum number of comparisons to perform
        
    Returns:
    --------
    list
        List of potential duplicate pairs
    """
    import time
    from collections import defaultdict
    
    logger.info(f"Finding potential duplicates with {method} method")
    logger.info(f"Similarity threshold: {similarity_threshold}, Clinical threshold: {clinical_threshold}")
    start_time = time.time()
    
    # Group patients by study
    patients_by_study = defaultdict(list)
    for patient_key, profile in patient_profiles.items():
        patients_by_study[profile['study_id']].append(patient_key)
    
    logger.info(f"Patients grouped into {len(patients_by_study)} studies")
    
    # Initialize list for potential duplicates
    potential_duplicates = []
    
    # Track comparison count
    comparison_count = 0
    
    # Compare patients across different studies
    studies = list(patients_by_study.keys())
    
    for i in range(len(studies)):
        for j in range(i+1, len(studies)):
            study1 = studies[i]
            study2 = studies[j]
            
            patients_study1 = patients_by_study[study1]
            patients_study2 = patients_by_study[study2]
            
            logger.info(f"Comparing {len(patients_study1)} patients from {study1} with "
                       f"{len(patients_study2)} patients from {study2}")
            
            # Limitation for large datasets
            if max_comparisons and comparison_count >= max_comparisons:
                logger.info(f"Reached maximum number of comparisons: {max_comparisons}")
                break
            
            # Compare patients between these two studies
            for patient1 in patients_study1:
                for patient2 in patients_study2:
                    comparison_count += 1
                    
                    # Get patient profiles
                    profile1 = patient_profiles[patient1]
                    profile2 = patient_profiles[patient2]
                    
                    details = {}
                    similarity = 0.0
                    
                    # Clinical similarity
                    if method in ['combined', 'clinical_only']:
                        clinical_score, clinical_details = compute_clinical_similarity(
                            profile1, profile2, min_matches=min_matches
                        )
                        
                        # For clinical_only method, use the clinical score
                        if method == 'clinical_only':
                            similarity = clinical_score
                            details = clinical_details
                    
                    # Mutation similarity
                    if method in ['combined', 'mutations_only'] and mutation_profiles:
                        # Check if both patients have mutation profiles
                        if patient1 in mutation_profiles and patient2 in mutation_profiles:
                            mut_profile1 = mutation_profiles[patient1]
                            mut_profile2 = mutation_profiles[patient2]
                            
                            mutation_score, mutation_details = compute_mutation_similarity(
                                mut_profile1, mut_profile2
                            )
                            
                            # For mutations_only method, use the mutation score
                            if method == 'mutations_only':
                                similarity = mutation_score
                                details = mutation_details
                        else:
                            # If mutation profiles missing, set defaults
                            if method == 'mutations_only':
                                similarity = 0.0
                                details = {'error': 'Missing mutation data'}
                            mutation_score = 0.0
                            mutation_details = {'error': 'Missing mutation data'}
                    else:
                        # If method is mutations_only but no mutation profiles provided
                        if method == 'mutations_only':
                            similarity = 0.0
                            details = {'error': 'No mutation data provided'}
                        mutation_score = 0.0
                        mutation_details = {'error': 'No mutation data provided'}
                    
                    # For combined method, calculate the combined score
                    if method == 'combined':
                        # Only proceed to combined if clinical score is high enough
                        if clinical_score >= clinical_threshold / 2:  # Lower threshold for combined comparison
                            similarity, details = compute_combined_similarity(
                                clinical_score, clinical_details,
                                mutation_score, mutation_details
                            )
                        else:
                            similarity = 0.0
                            details = {'error': 'Clinical score below threshold for combined evaluation'}
                    
                    # Different thresholds based on method
                    threshold = similarity_threshold
                    if method == 'clinical_only':
                        threshold = clinical_threshold
                    
                    # If similarity exceeds threshold, consider as potential duplicate
                    if similarity >= threshold:
                        # Get clinical matching details if available
                        clinical_info = ""
                        if isinstance(details.get('clinical', {}), dict) and 'field_details' in details['clinical']:
                            matches = []
                            for field, match_score in details['clinical']['field_details'].items():
                                if match_score > 0.5:  # Consider substantial matches
                                    matches.append(field)
                            if matches:
                                clinical_info = ", ".join(matches)
                        
                        potential_duplicates.append({
                            'patient1': patient1,
                            'patient2': patient2,
                            'study1': study1,
                            'study2': study2,
                            'similarity': similarity,
                            'details': details,
                            'matching_fields': clinical_info
                        })
    
    logger.info(f"Completed {comparison_count} comparisons in {time.time() - start_time:.2f} seconds")
    logger.info(f"Found {len(potential_duplicates)} potential duplicate pairs")
    
    # Sort by similarity (highest first)
    potential_duplicates.sort(key=lambda x: x['similarity'], reverse=True)
    
    return potential_duplicates