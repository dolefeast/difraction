import matplotlib.pyplot as plt
import scipy.fft as fft
import numpy as np


class Monochromatic_Experiment:
    def __init__(self, wavelength, N, size):
        self.N = N #Number of pixels in the slit
        self.size = size #size of the slit, in mm
        self.wavelength = wavelength #Wavelength of the incoming wave

        #Coordinate space coordinates
        self.x = np.linspace(0, size, N)
        self.y = np.linspace(0, size, N)
        self.X, self.Y = np.meshgrid(self.x, self.y)

        #Frequency space coordinates
        kx = np.linspace(-2*np.pi * N//2/size, 2*np.pi*N//2 /size, N)
        ky = np.linspace(-2*np.pi * N//2 /size, 2*np.pi*N//2 /size, N)
        self.Kx, self.Ky = np.meshgrid(kx, ky)

        self.slit = np.zeros((N, N)) + np.zeros((N, N))*1j #So we can add complex functions

    def spherical_wave(self, r0):
        """In case the wave comes from a source at a distance r0<\inf
        """
        E = np.exp(1j * 2 * np.pi/self.wavelength * np.sqrt(self.X**2 + self.Y**2+ r0**2))
        self.slit *= E

    def rectangle(self, width, height, x0, y0):
        """Adds a rectangular slit to the setup. All units in mm.
        Params:
        width: int
            Width in mm of the rectangular slit
        height: int
            Height in mm of the rectangular slit
        x0: int
            Horizontal position of the center of the recangular slit
        y0: int
            Vertical position of the center of the recangular slit
        """
        #second try: scaling
        self.scale_factor = self.N//self.size
        x0 *= self.scale_factor
        y0 *= self.scale_factor
        height *= self.scale_factor
        width *= self.scale_factor
        self.slit[x0 - height//2:x0 + height//2, int(y0 - width//2):int(y0 + width//2)] = 1+0j 
        
        #self.slit[x0 - height//2:x0 + height//2, y0 - width//2:y0 + width//2] = 1+0j 

    def grid(self, width, height, x0, y0, N, offset):
        for n in range(N):
            exp.rectangle(width, height, self.N//2, self.N//2 - n*offset)
            exp.rectangle(width, height, self.N//2, self.N//2 + n*offset)
            exp.rectangle(height, width, self.N//2 - n*offset, self.N//2 )
            exp.rectangle(height, width, self.N//2 + n*offset, self.N//2 )

    def propagate(self, Z, A):
        """
        Params:
        Z: The distance from the screen to the slit
        A: The Fourier transform A(kx, ky, 0)

        Returns the forier transform A(kx, ky, Z)
        """
        phase = np.exp(1j*Z*np.sqrt( (2*np.pi/self.wavelength)**2 - self.Kx**2 - self.Ky**2)) 
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

N = 200 # pixels number
size = 10 * mm # physical length of a square (size_x = size_y)

wavelength = 900*nm

offset = 0.2 #offset of the slits. Half the distance between the centers of the slits.
exp = Monochromatic_Experiment(wavelength, N, size, r0)

#exp.grid(2, 100, N//2, N//2, 30, 20)
exp.rectangle(0.1, 9, size//2, size//2 + offset)
exp.rectangle(0.1, 9, size//2, size//2 - offset)
#exp.rectangle(50, 50, N//2, N//2)
exp.spherical_wave(20)

fig, ax = plt.subplots()

ax.tick_params(left=False,
                bottom=True,
                labelleft=False,
                labelbottom=True)

ax.imshow(exp.field_at_Z(1000))

for field, Z in exp.evolve(0, 1000, 300):  
    ax.cla()
    ax.set_title('Distance: ' +str(int(Z)) + ' mm\nSide= ' + str(int(size)) + ' mm')
    exp.r0 += Z 
    ax.imshow(field)
    plt.pause(0.0000001)
plt.show()
