import os
import numpy as np
from dxtbx.model import ExperimentList

class ExptList():
    """
    Reads a DIALS experiment file and contains functions for parsing data inside.
    """

    #-------------------------------------------------------------------#
    # Constructor
    
    def __init__(self, exptfile):
        self.readExptFile(exptfile)

    #-------------------------------------------------------------------#
    # Attributes
    
    @property
    def elist(self):
        return self._elist

    @property
    def crystals(self):
        return self._crystals

    @property
    def images(self):
        return self._images

    #-------------------------------------------------------------------#
    # I/O Methods
    
    def readExptFile(self, exptfile):
        """
        Read Precognition .inp file and update geometric attributes

        Parameters
        ----------
        exptfile : str
            Path to .expt file from which to read.

        Notes
        -----
        It is assumed that the experiments are stills.
        """
        # Check that inpfile exists
        if not os.path.exists(exptfile):
            raise ValueError(f"Cannot find file: {exptfile}")

        # Create ExperimentList from file
        elist = ExperimentList.from_file(exptfile, check_format=False)

        # Check that experiments are stills
        assert len(elist.crystals()) == len(elist), 'ERROR: DIALS experiment file does not have stills. Please only use stills.'

        # Set attributes
        self._elist = elist
        self._crystals = elist.crystals()
        self._images = elist.imagesets()
        return
            
    #-------------------------------------------------------------------#
    # Crystallographic Methods
    def get_image_filenames(self):
        imgs = self.images
        image_filenames = []
        for img in imgs:
            image_filenames.append(img.paths()[0])
        return image_filenames

    def get_orthogonalization_matrix(self):
        """
        Compute real-space orthogonalization matrix from unit cell parameters.
        """
        cryst = self.crystals[0]
        cell = cryst.get_unit_cell()
        O = np.reshape(cell.orthogonalization_matrix(), (3,3))
        return O.T

    def get_reciprocal_Amatrices(self):
        """
        Get A matrix in reciprocal lattice basis (A*)
        """
        crystals = self.crystals
        A_stars = []
        for cryst in crystals:
            A_stars.append(np.reshape(cryst.get_A(), (3,3)))
        return A_stars
