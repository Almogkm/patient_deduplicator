"""
Input/Output utilities for the Patient Deduplicator package.
"""

import time
import logging
import pandas as pd
from pathlib import Path

logger = logging.getLogger('patient_deduplication.io')

def load_data(mutations_file, patients_file, sample_size=None):
    """
    Load the mutations and patients data from parquet files.
    
    Parameters:
    -----------
    mutations_file : str
        Path to the mutations data file (parquet format)
    patients_file : str
        Path to the patients clinical data file (parquet format)
    sample_size : int, optional
        If provided, only load a random sample of this size
    
    Returns:
    --------
    tuple
        (mutations_dataframe, patients_dataframe)
    """
    start_time = time.time()
    
    # Load mutations data
    logger.info(f"Loading mutations data from {mutations_file}")
    if sample_size:
        # When working with large files, sometimes sampling is useful for initial exploration
        mutations_df = pd.read_parquet(mutations_file).sample(sample_size)
    else:
        mutations_df = pd.read_parquet(mutations_file)
    
    # Log basic statistics about the mutations dataset
    logger.info(f"Mutations data loaded in {time.time() - start_time:.2f} seconds")
    logger.info(f"Mutations dataset shape: {mutations_df.shape}")
    
    # Check if required columns are present in mutations data
    mutations_required_cols = ['hugoGeneSymbol', 'uniquePatientKey', 'patientId', 'studyId', 
                     'mutationType', 'proteinChange', 'chr', 'startPosition']
    
    mutations_missing_cols = [col for col in mutations_required_cols if col not in mutations_df.columns]
    if mutations_missing_cols:
        logger.warning(f"Missing required columns in mutations data: {mutations_missing_cols}")
    
    # Load patients data
    patients_start_time = time.time()
    logger.info(f"Loading patients data from {patients_file}")
    patients_df = pd.read_parquet(patients_file)
    
    # Log basic statistics about the patients dataset
    logger.info(f"Patients data loaded in {time.time() - patients_start_time:.2f} seconds")
    logger.info(f"Patients dataset shape: {patients_df.shape}")
    
    # Check if patient identifiers are present in patients data
    patients_required_cols = ['studyId', 'uniquePatientKey', 'patientId', 'SAMPLE_COUNT', 
                             'SEX', 'OS_STATUS', 'OS_MONTHS', 'RACE', 'ETHNICITY', 'AGE', 'GENDER']
    patients_missing_cols = [col for col in patients_required_cols if col not in patients_df.columns]
    if patients_missing_cols:
        logger.warning(f"Missing required columns in patients data: {patients_missing_cols}")
        
    # Check for overlapping patients between mutations and clinical data
    mutations_patients = set(mutations_df['uniquePatientKey'].unique())
    clinical_patients = set(patients_df['uniquePatientKey'].unique())
    overlap_count = len(mutations_patients.intersection(clinical_patients))
    
    logger.info(f"Found {overlap_count} patients with both mutation and clinical data")
    
    # Subset both dataframes to only include matched patients
    matched_patients = mutations_patients.intersection(clinical_patients)
    
    if matched_patients:
        # Subset mutations data
        original_mutations_count = len(mutations_df)
        mutations_df = mutations_df[mutations_df['uniquePatientKey'].isin(matched_patients)]
        mutations_removed = original_mutations_count - len(mutations_df)
        logger.info(f"Removed {mutations_removed} mutations for patients without clinical data")
        
        # Subset patients data
        original_patients_count = len(patients_df)
        patients_df = patients_df[patients_df['uniquePatientKey'].isin(matched_patients)]
        patients_removed = original_patients_count - len(patients_df)
        logger.info(f"Removed {patients_removed} patients without mutation data")
        
        # Log final dataset sizes
        logger.info(f"Final mutations dataset shape: {mutations_df.shape}")
        logger.info(f"Final patients dataset shape: {patients_df.shape}")
    else:
        logger.warning("No matched patients found between mutations and clinical data")
    
    return mutations_df, patients_df

def save_results(potential_duplicates, patient_profiles, output_path):
    """
    Export potential duplicate patients to a CSV file, including clinical match information.
    
    Parameters:
    -----------
    potential_duplicates : list
        List of dictionaries containing potential duplicate pairs
    patient_profiles : dict
        Dictionary of patient profiles
    output_path : str or Path
        Path to save the output file
    
    Returns:
    --------
    str
        Path to the saved file
    """
    if not potential_duplicates:
        logger.warning("No potential duplicates to export")
        return None
    
    # Create a dataframe from the potential duplicates, including clinical data
    export_records = []
    
    for dup in potential_duplicates:
        # Extract profile data for both patients
        profile1 = patient_profiles[dup['patient1']]
        profile2 = patient_profiles[dup['patient2']]
        
        # Basic record with identifiers and similarity
        record = {
            'patient1': dup['patient1'],
            'patient1_id': profile1.get('patientId', ''),
            'study1': dup['study1'],
            'patient2': dup['patient2'],
            'patient2_id': profile2.get('patientId', ''),
            'study2': dup['study2'],
            'similarity_score': dup['similarity'],
            'matching_fields': dup.get('matching_fields', '')
        }
        
        # Add mutation similarity details if available
        if 'mutation' in dup['details']:
            mutation_details = dup['details']['mutation']
            if isinstance(mutation_details, dict) and 'error' not in mutation_details:
                record.update({
                    'common_mutations': mutation_details.get('common_mutations', 0),
                    'total_mutations': mutation_details.get('total_unique_mutations', 0),
                    'common_genes': mutation_details.get('common_genes', 0),
                    'mutation_similarity': mutation_details.get('mutation_signature_similarity', 0)
                })
        
        # Add demographic comparison
        demographic_fields = ['SAMPLE_COUNT', 'SEX', 'OS_STATUS', 'OS_MONTHS', 'RACE', 'ETHNICITY', 'AGE', 'GENDER']
        
        for field in demographic_fields:
            field_lower = field.lower()
            
            # Check both upper and lower case field names
            value1 = profile1.get(field, profile1.get(field_lower, 'NA'))
            value2 = profile2.get(field, profile2.get(field_lower, 'NA'))
            
            record[f'{field}_patient1'] = value1
            record[f'{field}_patient2'] = value2
        
        export_records.append(record)
    
    # Create DataFrame
    df_export = pd.DataFrame(export_records)
    
    # Save to CSV
    output_path = Path(output_path)
    df_export.to_csv(output_path, index=False)
    logger.info(f"Exported potential duplicates to {output_path}")
    
    return str(output_path)