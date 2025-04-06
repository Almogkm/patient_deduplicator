"""
Preprocessing utilities for the Patient Deduplicator package.
"""

import logging
import pandas as pd
import numpy as np

logger = logging.getLogger('patient_deduplication.preprocessing')

def preprocess_patient_data(patients_df):
    """
    Preprocess the patients data for deduplication analysis.
    
    This includes:
    - Handling missing values
    - Preprocessing clinical variables for comparison
    - Creating additional features for comparison
    
    Parameters:
    -----------
    patients_df : pandas.DataFrame
        DataFrame containing patient clinical data
    
    Returns:
    --------
    pandas.DataFrame
        Preprocessed patient data
    """
    if patients_df is None or patients_df.empty:
        logger.error("No patient data to preprocess")
        return None
    
    logger.info("Preprocessing patient data")
    logger.info(f"Original patient data shape: {patients_df.shape}")
    
    # Make a copy to avoid modifying the original
    preprocessed_df = patients_df.copy()
    
    # Identify important clinical columns for patient matching
    demographics_cols = []
    for col in preprocessed_df.columns:
        # Add key demographic columns if they exist
        if col.upper() in ['AGE', 'SEX', 'GENDER', 'RACE', 'ETHNICITY']:
            demographics_cols.append(col)
    
    logger.info(f"Identified {len(demographics_cols)} demographic columns for matching: {demographics_cols}")
    
    # Normalize text fields (lowercase, strip whitespace)
    for col in preprocessed_df.columns:
        if preprocessed_df[col].dtype == 'object':
            preprocessed_df[col] = preprocessed_df[col].astype(str).str.strip().str.lower()
            
    # Handle common missing value patterns
    na_variants = ['unknown', 'unkown', 'u', 'n/a', 'na', 'nan', '']
    for variant in na_variants:
        for col in preprocessed_df.columns:
            if preprocessed_df[col].dtype == 'object':
                preprocessed_df.loc[preprocessed_df[col] == variant, col] = np.nan
    
    # For age columns, ensure numeric type
    age_cols = [col for col in preprocessed_df.columns if 'age' in col.lower()]
    for col in age_cols:
        if col in preprocessed_df.columns:
            try:
                preprocessed_df[col] = pd.to_numeric(preprocessed_df[col], errors='coerce')
            except:
                logger.warning(f"Could not convert {col} to numeric")
    
    # For categorical demographics, standardize common values
    if 'SEX' in preprocessed_df.columns or 'GENDER' in preprocessed_df.columns:
        # Standardize sex/gender
        sex_col = 'SEX' if 'SEX' in preprocessed_df.columns else 'GENDER'
        sex_map = {
            'male': 'male', 'm': 'male',
            'female': 'female', 'f': 'female',
            'other': 'other'
        }
        preprocessed_df[sex_col] = preprocessed_df[sex_col].map(lambda x: sex_map.get(str(x).lower(), x))
    
    logger.info(f"Preprocessed patient data shape: {preprocessed_df.shape}")
    return preprocessed_df

def create_patient_profiles(patients_df):
    """
    Create basic patient profiles based on clinical data.
    
    Parameters:
    -----------
    patients_df : pandas.DataFrame
        DataFrame containing preprocessed patient clinical data
    
    Returns:
    --------
    dict
        Dictionary of patient profiles keyed by uniquePatientKey
    """
    if patients_df is None or patients_df.empty:
        logger.error("No preprocessed data available for creating profiles")
        return {}
    
    logger.info("Creating basic patient profiles")
    
    # Initialize patient profiles dictionary
    patient_profiles = {}
    
    # Define demographic columns to extract
    demographic_columns = []
    for col in patients_df.columns:
        # Standard demographic columns
        if col.upper() in ['AGE', 'SEX', 'GENDER', 'RACE', 'ETHNICITY']:
            demographic_columns.append(col)
    
    # Process clinical data
    clinical_patients = patients_df['uniquePatientKey'].unique()
    logger.info(f"Creating profiles for {len(clinical_patients)} patients")
    
    # Track progress
    processed_count = 0
    report_interval = max(1, len(clinical_patients) // 10)  # Report ~10 times
    
    for patient_key in clinical_patients:
        try:
            # Get patient clinical data
            patient_clinical = patients_df[patients_df['uniquePatientKey'] == patient_key]
            
            if patient_clinical.empty:
                continue
            
            # Create basic profile
            patient_profiles[patient_key] = {
                'patient_id': patient_key,
                'study_id': patient_clinical['studyId'].iloc[0],
                'has_clinical_data': True,
                'has_mutation_data': False  # Will be updated later for similar patients
            }
            
            # Add demographic and clinical information
            for col in demographic_columns:
                if col in patient_clinical.columns and not patient_clinical[col].isna().all():
                    # Extract the value, handling both single and multiple rows
                    if len(patient_clinical) == 1:
                        value = patient_clinical[col].iloc[0]
                    else:
                        # If multiple rows, get most common non-null value
                        value_counts = patient_clinical[col].value_counts(dropna=True)
                        if not value_counts.empty:
                            value = value_counts.index[0]
                        else:
                            value = None
                            
                    # Add to profile if value exists
                    if value is not None:
                        patient_profiles[patient_key][col] = value
            
            # Update progress
            processed_count += 1
            if processed_count % report_interval == 0:
                logger.info(f"Processed {processed_count}/{len(clinical_patients)} patients")
                
        except Exception as e:
            logger.error(f"Error processing clinical data for patient {patient_key}: {str(e)}")
    
    logger.info(f"Created basic profiles for {len(patient_profiles)} patients")
    return patient_profiles

def create_mutation_profiles(mutations_df, patient_keys=None):
    """
    Create mutation profiles for the specified patients.
    
    Parameters:
    -----------
    mutations_df : pandas.DataFrame
        DataFrame containing mutation data
    patient_keys : list, optional
        List of patient keys to create profiles for. If None, create for all patients.
    
    Returns:
    --------
    dict
        Dictionary mapping patient keys to their mutation profiles
    """
    if mutations_df is None or mutations_df.empty:
        logger.error("No mutation data available")
        return {}
        
    # Filter mutations for these patients if specified
    if patient_keys:
        filtered_mutations = mutations_df[mutations_df['uniquePatientKey'].isin(patient_keys)]
        logger.info(f"Creating mutation profiles for {len(patient_keys)} specified patients")
    else:
        filtered_mutations = mutations_df
        patient_keys = filtered_mutations['uniquePatientKey'].unique()
        logger.info(f"Creating mutation profiles for all {len(patient_keys)} patients")
    
    # Group by patient
    mutation_groups = filtered_mutations.groupby('uniquePatientKey')
    
    # Create mutation profiles
    mutation_profiles = {}
    for patient_key, mutations in mutation_groups:
        # Create mutation signature (gene + protein change + chr + position)
        if 'proteinChange' in mutations.columns and 'chr' in mutations.columns and 'startPosition' in mutations.columns:
            mutations['mutation_signature'] = mutations.apply(
                lambda row: f"{row['hugoGeneSymbol']}:{row.get('proteinChange', '')}:{row['chr']}:{row['startPosition']}", 
                axis=1
            )
        else:
            # Fallback if not all columns are available
            mutations['mutation_signature'] = mutations['hugoGeneSymbol']
            
        # Build profile
        mutation_profiles[patient_key] = {
            'mutation_signatures': set(mutations['mutation_signature']),
            'mutated_genes': set(mutations['hugoGeneSymbol']),
            'total_mutations': len(mutations),
            'mutation_types': mutations['mutationType'].value_counts().to_dict() if 'mutationType' in mutations.columns else {},
            'chr_distribution': mutations['chr'].value_counts().to_dict() if 'chr' in mutations.columns else {}
        }
            
    logger.info(f"Created mutation profiles for {len(mutation_profiles)} patients")
    return mutation_profiles