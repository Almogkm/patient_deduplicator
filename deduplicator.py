"""
Patient Deduplicator - Core module

This module contains the main PatientDeduplicator class responsible for identifying
potential duplicate patients across different studies based on clinical and genetic data.
"""

import time
import logging
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from .utils.preprocessing import preprocess_patient_data, create_patient_profiles, create_mutation_profiles
from .utils.similarity import compute_clinical_similarity, compute_mutation_similarity, find_potential_duplicates
from .utils.io import load_data, save_results
from .utils.visualization import plot_potential_duplicates, plot_similarity_distribution, plot_study_heatmap

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('patient_deduplication')


class PatientDeduplicator:
    """
    A class to identify potential duplicate patients across different studies
    based on both genetic mutation profiles and clinical patient characteristics.
    """
    
    def __init__(self, mutations_file=None, patients_file=None, output_dir=None):
        """
        Initialize the deduplicator with mutations and patients file paths.
        
        Parameters:
        -----------
        mutations_file : str, optional
            Path to the mutations data file (parquet format)
        patients_file : str, optional
            Path to the patients clinical data file (parquet format)
        output_dir : str, optional
            Directory to save output files
        """
        self.mutations_file = mutations_file
        self.patients_file = patients_file
        self.output_dir = Path(output_dir) if output_dir else Path.cwd() / "output"
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        # Data containers
        self.mutations_df = None
        self.patients_df = None
        self.patient_profiles = {}
        self.mutation_profiles = {}
        self.potential_duplicates = []
        
        # Configuration settings
        self.similarity_threshold = 0.7
        self.clinical_threshold = 0.8
        self.min_matches = 2
        self.method = 'combined'
        
        logger.info(f"Initialized PatientDeduplicator")
        if mutations_file and patients_file:
            logger.info(f"- Mutations file: {mutations_file}")
            logger.info(f"- Patients file: {patients_file}")
        logger.info(f"Output will be saved to: {self.output_dir}")
    
    def set_files(self, mutations_file, patients_file):
        """
        Set or update the file paths.
        
        Parameters:
        -----------
        mutations_file : str
            Path to the mutations data file (parquet format)
        patients_file : str
            Path to the patients clinical data file (parquet format)
        """
        self.mutations_file = mutations_file
        self.patients_file = patients_file
        logger.info(f"Updated file paths:")
        logger.info(f"- Mutations file: {mutations_file}")
        logger.info(f"- Patients file: {patients_file}")
        
    def set_configuration(self, similarity_threshold=None, clinical_threshold=None, 
                        min_matches=None, method=None):
        """
        Configure the deduplication parameters.
        
        Parameters:
        -----------
        similarity_threshold : float, optional
            Threshold for considering patients as potential duplicates (0-1)
        clinical_threshold : float, optional
            Higher threshold for clinical-only matches (0-1)
        min_matches : int, optional
            Minimum number of matching fields required for a valid comparison
        method : str, optional
            Similarity method to use ('combined', 'mutations_only', 'clinical_only')
        """
        if similarity_threshold is not None:
            self.similarity_threshold = similarity_threshold
        
        if clinical_threshold is not None:
            self.clinical_threshold = clinical_threshold
            
        if min_matches is not None:
            self.min_matches = min_matches
            
        if method is not None:
            if method in ['combined', 'mutations_only', 'clinical_only']:
                self.method = method
            else:
                logger.warning(f"Invalid method: {method}. Using '{self.method}' instead.")
        
        logger.info(f"Configuration updated:")
        logger.info(f"- Similarity threshold: {self.similarity_threshold}")
        logger.info(f"- Clinical threshold: {self.clinical_threshold}")
        logger.info(f"- Minimum matches: {self.min_matches}")
        logger.info(f"- Method: {self.method}")
        
    def load_data(self, sample_size=None):
        """
        Load the mutations and patients data from parquet files.
        
        Parameters:
        -----------
        sample_size : int, optional
            If provided, only load a random sample of this size
        
        Returns:
        --------
        tuple
            (mutations_dataframe, patients_dataframe)
        """
        if not self.mutations_file or not self.patients_file:
            logger.error("Mutations or patients file path not set. Use set_files() method first.")
            return None, None
            
        self.mutations_df, self.patients_df = load_data(
            self.mutations_file, 
            self.patients_file, 
            sample_size=sample_size
        )
        
        return self.mutations_df, self.patients_df
    
    def preprocess_patient_data(self):
        """
        Preprocess the patients data for deduplication analysis.
        
        Returns:
        --------
        pandas.DataFrame
            Preprocessed patient data
        """
        if self.patients_df is None:
            logger.error("No patient data loaded. Call load_data() first.")
            return None
            
        logger.info("Preprocessing patient data")
        self.patients_df = preprocess_patient_data(self.patients_df)
        
        return self.patients_df
        
    def create_patient_profiles(self):
        """
        Create basic patient profiles based on clinical data only.
        
        Returns:
        --------
        dict
            Dictionary of patient profiles keyed by uniquePatientKey
        """
        if self.patients_df is None:
            logger.error("No preprocessed data available. Call preprocess_patient_data() first.")
            return None
            
        logger.info("Creating patient profiles")
        self.patient_profiles = create_patient_profiles(self.patients_df)
        
        return self.patient_profiles
        
    def create_mutation_profiles(self, patient_keys=None):
        """
        Create mutation profiles for patients.
        
        Parameters:
        -----------
        patient_keys : list, optional
            List of patient keys to create mutation profiles for. 
            If None, create profiles for all patients.
            
        Returns:
        --------
        dict
            Dictionary of mutation profiles keyed by uniquePatientKey
        """
        if self.mutations_df is None:
            logger.error("No mutation data loaded. Call load_data() first.")
            return None
            
        logger.info("Creating mutation profiles")
        self.mutation_profiles = create_mutation_profiles(
            self.mutations_df, 
            patient_keys=patient_keys
        )
        
        return self.mutation_profiles
        
    def find_potential_duplicates(self, max_comparisons=None):
        """
        Find potential duplicate patients across different studies.
        
        Parameters:
        -----------
        max_comparisons : int, optional
            Maximum number of comparisons to perform
            
        Returns:
        --------
        list
            List of potential duplicate pairs
        """
        if not self.patient_profiles:
            logger.error("No patient profiles available. Call create_patient_profiles() first.")
            return []
            
        logger.info(f"Finding potential duplicates with {self.method} method")
        
        self.potential_duplicates = find_potential_duplicates(
            self.patient_profiles,
            self.mutation_profiles,
            similarity_threshold=self.similarity_threshold,
            clinical_threshold=self.clinical_threshold,
            min_matches=self.min_matches,
            method=self.method,
            max_comparisons=max_comparisons
        )
        
        return self.potential_duplicates
        
    def export_results(self, filename=None):
        """
        Export potential duplicate patients to a CSV file.
        
        Parameters:
        -----------
        filename : str, optional
            Name of the output file. If None, a timestamped filename will be used.
            
        Returns:
        --------
        str
            Path to the saved file
        """
        if not self.potential_duplicates:
            logger.warning("No potential duplicates to export")
            return None
            
        # Generate filename if not provided
        if filename is None:
            timestamp = time.strftime("%Y%m%d-%H%M%S")
            filename = f"potential_duplicates_{timestamp}.csv"
            
        # Create full path
        output_path = self.output_dir / filename
        
        # Export results
        csv_path = save_results(
            self.potential_duplicates,
            self.patient_profiles,
            output_path
        )
        
        return csv_path
        
    def visualize_results(self, max_pairs=20):
        """
        Visualize the potential duplicate patient pairs.
        
        Parameters:
        -----------
        max_pairs : int
            Maximum number of pairs to visualize
            
        Returns:
        --------
        dict
            Dictionary with paths to visualization files
        """
        if not self.potential_duplicates:
            logger.warning("No potential duplicates to visualize")
            return None
            
        logger.info("Generating visualizations")
        
        # Create visualizations
        viz_paths = {}
        
        # Main duplicate pairs visualization
        viz_paths['duplicates'] = plot_potential_duplicates(
            self.potential_duplicates,
            self.patient_profiles,
            self.output_dir,
            max_pairs=max_pairs
        )
        
        # Similarity distribution
        viz_paths['similarity_distribution'] = plot_similarity_distribution(
            self.potential_duplicates,
            self.output_dir
        )
        
        # Study heatmap
        viz_paths['study_heatmap'] = plot_study_heatmap(
            self.potential_duplicates,
            self.output_dir
        )
        
        return viz_paths
        
    def run_pipeline(self, sample_size=None, max_comparisons=None, visualize=True):
        """
        Run the complete patient deduplication pipeline.
        
        Parameters:
        -----------
        sample_size : int, optional
            If provided, only load a random sample of this size
        max_comparisons : int, optional
            Maximum number of comparisons to perform
        visualize : bool
            Whether to generate visualizations
            
        Returns:
        --------
        dict
            Results dictionary with potential duplicates and output paths
        """
        start_time = time.time()
        logger.info("Starting patient deduplication pipeline")
        logger.info(f"Configuration: method={self.method}, min_matches={self.min_matches}")
        
        # Step 1: Load data
        self.load_data(sample_size=sample_size)
        
        # Step 2: Preprocess patient data
        self.preprocess_patient_data()
        
        # Step 3: Create patient profiles
        self.create_patient_profiles()
        
        # Step 4: Create mutation profiles (if needed)
        if self.method in ['combined', 'mutations_only']:
            self.create_mutation_profiles()
        
        # Step 5: Find potential duplicates
        self.find_potential_duplicates(max_comparisons=max_comparisons)
        
        # Step 6: Export results
        csv_path = self.export_results()
        
        # Step 7: Visualize results (if requested)
        viz_paths = None
        if visualize:
            viz_paths = self.visualize_results()
        
        logger.info(f"Pipeline completed in {time.time() - start_time:.2f} seconds")
        
        return {
            'potential_duplicates': self.potential_duplicates,
            'csv_output': csv_path,
            'visualizations': viz_paths
        }