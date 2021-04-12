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
        self.X, Y = np.meshgrid(self.x, self.y)

        #Frequency space coordinates
        kx = np.linspace(-2*np.pi * N, -2*np.pi, N)
        ky = np.linspace(-2*np.pi * N, -2*np.pi, N)
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

        fourier_at_0 = self.fourier_transform()
        fourier_at_z = fft.fftshift(self.propagate(Z, fourier_at_0))
        field = fft.ifft2(fourier_at_z)

        return self.intensity(field)





#Units definition
mm = 1
m = 1000*mm
cm = 10*m
nm = 1e-6 * mm

N = 200 # pixels number
size = 500 * mm # physical length of a square (extent_x = extent_y)


wavelength = 600*nm
Z = 1*m

exp = Monochromatic_Experiment(wavelength, N, size)
exp.rectangle(30, 30, 100, 100)




fig, ax = plt.subplots(1, 2)

ax[0].imshow(exp.slit)
for Z in range(100):
    print(Z)
    ax[1].cla()
    ax[1].imshow(exp.field_at_Z(50))
    plt.pause(0.01)
#plt.savefig('img.pdf')
plt.show()


