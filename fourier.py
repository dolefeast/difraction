import matplotlib.pyplot as plt
import scipy.fft as fft
import numpy as np


class Monochromatic_Experiment:
    def __init__(self, wavelength, N, size):
        self.N = N #Number of pixels in the slit
        self.size = size #size of the slit, in mm
        self.wavelength = wavelength #Wavelength of the incoming wave

        #Coordinate space coordinates
        self.x = np.linspace(-size/2, size/2, N)
        self.y = np.linspace(-size/2, size/2, N)
        self.X, self.Y = np.meshgrid(self.x, self.y)

        #Frequency space coordinates
        kx = np.linspace(-2*np.pi * N//2/size, 2*np.pi*N//2 /size, N)
        ky = np.linspace(-2*np.pi * N//2 /size, 2*np.pi*N//2 /size, N)
        self.Kx, self.Ky = np.meshgrid(kx, ky)

        self.slit = np.zeros((N, N))

    def rectangle(self, width, height, x0, y0):
        """Adds a rectangular slit to the setup. All units in px.
        Params:
        width: int
            Width in pixels of the rectangular slit
        height: int
            Height in pixels of the rectangular slit
        x0: int
            Horizontal position of the center of the recangular slit
        y0: int
            Vertical position of the center of the recangular slit
        """
        self.slit[x0 - height//2:x0 + height//2, y0 - width//2:y0 + width//2] += 1

    def propagate(self, Z, A):
        """
        Params:
        Z: The distance from the screen to the slit
        A: The Fourier transform A(kx, ky, 0)

        Returns the forier transform A(kx, ky, Z)
        """
        phase = np.exp(1j*Z*np.sqrt( (2*np.pi/self.wavelength)**2 - self.Kx**2 - self.Ky**2)) #Not 100% sure of this
        return A*phase

    def fourier_transform(self):
        """Returns the fourier transform of the slit.
        """
        fourier = fft.fft2(self.slit)
        return fft.fftshift(fourier)

    def intensity(self, F):
        """Returns the intensity of a field F
        Params:
        F: float
        """
        return np.real(F*np.conjugate(F))

    def field_at_Z(self, Z):
        """Returns the interference pattern observed at a distance Z from the slit. Kind of a 'full experiment'
        Params:
        Z: float
        Distance, in mm, from the slit to the 'screen'
        Returns:
        self.N x self.N array of intensities to be plotted"""

        fourier_at_0 = self.fourier_transform() #Calculate the Fourier transform of the signal (Field at z=0)
        fourier_at_Z = self.propagate(Z, fourier_at_0) #The Fourier transform of the field at z = Z.
        #fourier_at_Z = fft.fftshift(fourier_at_Z) #Undo the shift

        field = fft.ifft2(fourier_at_Z)
        
        return self.intensity(field)
    def evolve(self, Zmin, Zmax, N=100):
        """Returns a list of the fields for different Z"""
        for Z in np.linspace(Zmin, Zmax, N): #A stupidly complex way of changing Z
            yield (self.field_at_Z(Z), Z)


#Units definition
mm = 1
m = 1000*mm
cm = 10*m
nm = 1e-6 * mm

N = 500 # pixels number
size = 20 * mm # physical length of a square (extent_x = extent_y)


wavelength = 600*nm

A = 10 #offset of the slits. Half the distance between the centers of the slits.
exp = Monochromatic_Experiment(wavelength, N, size)
exp.rectangle(4, 100, N//2, N//2 - A)
exp.rectangle(4, 100, N//2, N//2 + A)


fig, ax = plt.subplots(1, 2)

ax[0].imshow(exp.slit)

for field, Z in exp.evolve(0, 400, 40):  
    ax[1].cla()
    ax[1].set_title(str(Z) + 'mm')
    ax[1].imshow(field)
    plt.pause(0.000001)
plt.show()


