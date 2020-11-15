import os
import numpy as np

class FrameGeometry():
    """
    Provides a representation of  per-frame geometry in Precognition.

    Includes attributes for accessing geometric parameters, and methods
    for reading and writing Precognition .inp geometry files. 
    """

    #-------------------------------------------------------------------#
    # Constructor
    
    def __init__(self, inpfile):
        self.goniometer = None
        self.image = None
        self.readINPFile(inpfile)

    #-------------------------------------------------------------------#
    # Attributes
        
    @property
    def crystal(self):
        return self._crystal

    @property
    def a(self):
        return self._crystal[0]

    @property
    def b(self):
        return self._crystal[1]

    @property
    def c(self):
        return self._crystal[2]

    @property
    def alpha(self):
        return self._crystal[3]

    @property
    def beta(self):
        return self._crystal[4]

    @property
    def gamma(self):
        return self._crystal[5]
    
    @crystal.setter
    def crystal(self, values):
        if len(values) != 6:
            raise ValueError(f"Cell parameters must have 6 values")
        values = *map(float, values),
        self._crystal = values
        return

    @property
    def spacegroup(self):
        return self._spacegroup

    @spacegroup.setter 
    def spacegroup(self, value):
        self._spacegroup = value
        return

    @property
    def matrix(self):
        return self._matrix

    @matrix.setter
    def matrix(self, values):
        self._matrix = values
        return

    @property
    def omega(self):
        return self._omega

    @omega.setter
    def omega(self, values):
        self._omega = values
        return

    @property
    def goniometer(self):
        return self._goniometer

    @goniometer.setter
    def goniometer(self, values):
        self._goniometer = values
        return

    @property
    def imageformat(self):
        return self._imageformat

    @imageformat.setter
    def imageformat(self, value):
        self._imageformat = value
        return

    @property
    def distance(self):
        return self._distance

    @distance.setter
    def distance(self, values):
        self._distance = values
        return

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, values):
        self._center = values
        return

    @property
    def pixel(self):
        return self._pixel

    @pixel.setter
    def pixel(self, values):
        self._pixel = values
        return

    @property
    def swing(self):
        return self._swing

    @swing.setter
    def swing(self, values):
        self._swing = values
        return

    @property
    def tilt(self):
        return self._tilt

    @tilt.setter
    def tilt(self, values):
        self._tilt = values
        return
    
    @property
    def bulge(self):
        return self._bulge

    @bulge.setter
    def bulge(self, values):
        self._bulge = values
        return

    @property
    def image(self):
        return self._image

    @image.setter
    def image(self, value):
        self._image = value
        return

    @property
    def resolution(self):
        return self._resolution

    @resolution.setter
    def resolution(self, values):
        self._resolution = values
        return

    @property
    def wavelength(self):
        return self._wavelength

    @wavelength.setter
    def wavelength(self, values):
        self._wavelength = values
        return

    #-------------------------------------------------------------------#
    # I/O Methods
    
    def readINPFile(self, inpfile):
        """
        Read Precognition .inp file and update geometric attributes

        Parameters
        ----------
        inpfile : str
            Path to .inp file from which to read

        Notes
        -----
        It is assumed that the format of the .inp file is as follows:

        Input
           Field1     Values1
           Field2     Values2
           ...
           Quit
        """
        # Check that inpfile exists
        if not os.path.exists(inpfile):
            raise ValueError(f"Cannot find file: {inpfile}")

        # Read inpfile
        with open(inpfile, "r") as inp:
            lines = inp.readlines()
        if not ("Input" in lines[0] and "Quit" in lines[-1]):
            raise ValueError(f"{inpfile} does not meet formatting assumptions")

        # Parse geometric parameters
        geometry = {}
        for l in lines[1:-1]:
            fields = l.split()
            if fields != []:
                geometry[fields[0]] = fields[1:]

        for key, values in geometry.items():
            if key == "Crystal":
                self.crystal = values[:-1]
                self.spacegroup = values[-1]
            elif key == "Matrix":
                self.matrix = values
            elif key == "Omega":
                self.omega = values
            elif key == "Goniometer":
                self.goniometer = values
            elif key == "Format":
                self.imageformat = values[0]
            elif key == "Distance":
                self.distance = values
            elif key == "Center":
                self.center = values
            elif key == "Pixel":
                self.pixel = values
            elif key == "Swing":
                self.swing = values
            elif key == "Tilt":
                self.tilt = values
            elif key == "Bulge":
                self.bulge = values
            elif key == "Image":
                self.image = values
            elif key == "Resolution":
                self.resolution = values
            elif key == "Wavelength":
                self.wavelength = values
            else:
                raise ValueError(f"Unexpected key {key} in {inpfile}")
                
        return
                
    def writeINPFile(self, inpfile):
        """
        Write Precognition .inp file containing experimental geometry

        Parameters
        ----------
        inpfile : str
            Path to .inp file to which to write
        """

        line_crystal    = " ".join(map(str, self.crystal))
        line_matrix     = " ".join(self.matrix)
        line_omega      = " ".join(self.omega)

        line_distance   = " ".join(self.distance)
        line_center     = " ".join(self.center)
        line_pixel      = " ".join(self.pixel)
        line_swing      = " ".join(self.swing)
        line_tilt       = " ".join(self.tilt)
        line_bulge      = " ".join(self.bulge)
        line_resolution = " ".join(self.resolution)
        line_wavelength = " ".join(self.wavelength)

        if self.goniometer:
            line_goniometer = " ".join(self.goniometer)
            
            inp = (f"Input\n"
                   f"   Crystal    {line_crystal} {self.spacegroup}\n"
                   f"   Matrix     {line_matrix}\n"
                   f"   Omega      {line_omega}\n"
                   f"   Goniometer {line_goniometer}\n\n"
                   f"   Format     {self.imageformat}\n"
                   f"   Distance   {line_distance}\n"
                   f"   Center     {line_center}\n"
                   f"   Pixel      {line_pixel}\n"
                   f"   Swing      {line_swing}\n"
                   f"   Tilt       {line_tilt}\n"
                   f"   Bulge      {line_bulge}\n\n"
                   f"   Image {self.image[0]}    {self.image[1]}\n"
                   f"   Resolution {line_resolution}\n"
                   f"   Wavelength {line_wavelength}\n"
                   f"   Quit\n"
            )

        else:

            inp = (f"Input\n"
                   f"   Crystal    {line_crystal} {self.spacegroup}\n"
                   f"   Matrix     {line_matrix}\n"
                   f"   Omega      {line_omega}\n"
                   f"   Format     {self.imageformat}\n"
                   f"   Distance   {line_distance}\n"
                   f"   Center     {line_center}\n"
                   f"   Pixel      {line_pixel}\n"
                   f"   Swing      {line_swing}\n"
                   f"   Tilt       {line_tilt}\n"
                   f"   Bulge      {line_bulge}\n\n"
                   f"   Resolution {line_resolution}\n"
                   f"   Wavelength {line_wavelength}\n"
                   f"   Quit\n"
            )

        with open(inpfile, "w") as outfile:
            outfile.write(inp)

        return
            
    #-------------------------------------------------------------------#
    # Crystallographic Methods

    def get_orthogonalization_matrix(self):
        """
        Compute real-space orthogonalization matrix from unit cell parameters.
        """
        a = float(self.a)
        b = float(self.b)
        c = float(self.c)
        alpha = np.deg2rad(float(self.alpha))
        beta = np.deg2rad(float(self.beta))
        gamma = np.deg2rad(float(self.gamma))

        # Compute unit cell volume
        V = (a*b*c*np.sqrt(1 - np.cos(alpha)**2 -
                           np.cos(beta)**2 -
                           np.cos(gamma)**2 +
                           2*np.cos(alpha)*np.cos(beta)*np.cos(gamma)))
    
        # Compute Cartesian orthogonalization matrix (Rupp, Page 746)
        O = np.zeros((3, 3))
        O[0, 0] = a
        O[0, 1] = b*np.cos(gamma)
        O[1, 1] = b*np.sin(gamma)
        O[0, 2] = c*np.cos(beta)
        O[1, 2] = c*(np.cos(alpha) - (np.cos(beta)*np.cos(gamma)))/np.sin(gamma)
        O[2, 2] = V/(a*b*np.sin(gamma))

        return O.T

    def get_missetting_matrix(self):
        """
        Get missetting matrix for FrameGeometry
        """
        return np.array(self.matrix, dtype=np.float).reshape(3, 3)

    def get_goniometer_rotation_matrix(self):
        """
        Get rotation matrix associated with goniometer settings of FrameGeometry
        """
        o1 = np.deg2rad(float(self.omega[0]))
        o2 = np.deg2rad(float(self.omega[1]))    
        gonio_phi = np.deg2rad(float(self.goniometer[2]))
        
        def get_rotation_matrix(axis, angle):
            u = axis
            sin,cos = np.sin(angle),np.cos(angle)
            return cos*np.eye(3) + sin*np.cross(u, -np.eye(3)) + (1. - cos)*np.outer(u, u)

        R = get_rotation_matrix(np.array([0., 0., -1.]),  o1)
        R = get_rotation_matrix(np.array([0., 1., 0.]), o2)@R
        R = get_rotation_matrix((R@np.array([0., 1., 0.])[:,None])[:,0], gonio_phi)@R
        return R

    def get_reciprocal_Amatrix(self):
        """
        Get A matrix in reciprocal lattice basis (A*)
        """
        O = self.get_orthogonalization_matrix()
        missetting = self.get_missetting_matrix()
        R = self.get_goniometer_rotation_matrix()
        precog2mosflm = np.array(
            [[  0,  0,  1],
             [  0, -1,  0],
             [  1,  0,  0]]
        )
        A_star = precog2mosflm@(R@missetting@np.linalg.inv(O))
        return A_star

    def get_realspace_Amatrix(self):
        """
        Get A matrix in realspace lattice basis (A)
        """
        return np.linalg.inv(self.get_reciprocal_Amatrix())

    def get_realspace_unitcell_vectors(self):
        """
        Get unit cell vectors in realspace lattice basis (a, b, c)
        """
        a, b, c = self.get_realspace_Amatrix()
        return a, b, c
