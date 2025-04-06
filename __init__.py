"""
Patient Deduplicator Package

A package for identifying potential duplicate patients across different studies
based on clinical characteristics and genetic mutation profiles.
"""

__version__ = '0.1.0'

from .deduplicator import PatientDeduplicator

__all__ = ['PatientDeduplicator']