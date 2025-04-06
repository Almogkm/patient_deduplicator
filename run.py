#!/usr/bin/env python3
"""
Patient Deduplicator - Command Line Runner

This script provides a command-line interface to run the patient deduplication
pipeline with various configuration options.
"""

import os
import sys
import argparse
import logging
from pathlib import Path

from patient_deduplicator import PatientDeduplicator

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Patient Deduplicator - Find potential duplicate patients across studies"
    )
    
    # Required file paths
    parser.add_argument(
        "--mutations", "-m", 
        required=True,
        help="Path to mutations data file (parquet format)"
    )
    parser.add_argument(
        "--patients", "-p", 
        required=True,
        help="Path to patients clinical data file (parquet format)"
    )
    parser.add_argument(
        "--output", "-o", 
        default="output",
        help="Directory to save output files (default: output)"
    )
    
    # Configuration options
    parser.add_argument(
        "--sample", "-s", 
        type=int, 
        default=None,
        help="Sample size (number of patients to load, leave empty for full dataset)"
    )
    parser.add_argument(
        "--similarity", 
        type=float, 
        default=0.7,
        help="Similarity threshold (0-1) for considering patients as potential duplicates (default: 0.7)"
    )
    parser.add_argument(
        "--clinical", 
        type=float, 
        default=0.8,
        help="Clinical threshold (0-1) for considering patients as potential duplicates (default: 0.8)"
    )
    parser.add_argument(
        "--min-matches", 
        type=int, 
        default=2,
        help="Minimum number of matching fields required for a valid comparison (default: 2)"
    )
    parser.add_argument(
        "--method", 
        choices=["combined", "mutations_only", "clinical_only"],
        default="combined",
        help="Similarity method to use (default: combined)"
    )
    parser.add_argument(
        "--max-comparisons", 
        type=int, 
        default=None,
        help="Maximum number of comparisons to perform (leave empty for no limit)"
    )
    parser.add_argument(
        "--no-visualization", 
        action="store_true",
        help="Disable visualization generation"
    )
    parser.add_argument(
        "--verbose", "-v", 
        action="store_true",
        help="Enable verbose logging"
    )
    
    return parser.parse_args()

def main():
    """Main execution function"""
    # Parse arguments
    args = get_args()
    
    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger('patient_deduplication')
    
    # Print configuration
    logger.info("Patient Deduplicator - Starting with configuration:")
    logger.info(f"- Mutations file: {args.mutations}")
    logger.info(f"- Patients file: {args.patients}")
    logger.info(f"- Output directory: {args.output}")
    logger.info(f"- Sample size: {'Full dataset' if args.sample is None else args.sample}")
    logger.info(f"- Similarity threshold: {args.similarity}")
    logger.info(f"- Clinical threshold: {args.clinical}")
    logger.info(f"- Minimum matches: {args.min_matches}")
    logger.info(f"- Method: {args.method}")
    logger.info(f"- Max comparisons: {'No limit' if args.max_comparisons is None else args.max_comparisons}")
    logger.info(f"- Visualization: {'Disabled' if args.no_visualization else 'Enabled'}")
    
    try:
        # Initialize deduplicator
        deduplicator = PatientDeduplicator(
            mutations_file=args.mutations,
            patients_file=args.patients,
            output_dir=args.output
        )
        
        # Set configuration
        deduplicator.set_configuration(
            similarity_threshold=args.similarity,
            clinical_threshold=args.clinical,
            min_matches=args.min_matches,
            method=args.method
        )
        
        # Run pipeline
        results = deduplicator.run_pipeline(
            sample_size=args.sample,
            max_comparisons=args.max_comparisons,
            visualize=not args.no_visualization
        )
        
        # Report results
        logger.info("\nPipeline completed successfully!")
        logger.info(f"Found {len(results['potential_duplicates'])} potential duplicate patients")
        logger.info(f"Results saved to {results['csv_output']}")
        
        if results['visualizations']:
            for viz_type, path in results['visualizations'].items():
                logger.info(f"{viz_type.capitalize()} visualization saved to {path}")
                
        return 0
        
    except Exception as e:
        logger.error(f"Error during execution: {str(e)}", exc_info=args.verbose)
        return 1

if __name__ == "__main__":
    sys.exit(main())