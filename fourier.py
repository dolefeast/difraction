import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import scipy.fft as fft
import numpy as np

#Units definition
mm = 1
m = 1000 * mm
cm = 10 * mm
nm = 1e-6 * mm

class Monochromatic_Experiment:
    def __init__(self, wavelength, N, size):
        self.N = N #Number of pixels in the slit
        self.size = size #size of the experiment, in mm
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
        d = np.sqrt(self.X**2 + self.Y**2+ r0**2) #The distance to the source of every point in the slit (matrix)
        E = np.exp(1j * 2 * np.pi/self.wavelength * d)/d #The electric field (complex) of every point in the slit (matrix)
        self.slit *= E 

    def rectangle(self, width, height, x0, y0):
        """Adds a rectangular slit to the setup. All units in mm.
        Params:
        width: int
            Width in mm of the rectangular slit
        height: int
            Height in mm of the rectangular slit
        x0: int
            Horizontal position of the center of the rectangular slit
        y0: int
            Vertical position of the center of the rectangular slit
        """
        #second try: scaling
        self.scale_factor = self.N//self.size #To convert from mm to px
        x0 *= self.scale_factor
        y0 *= self.scale_factor
        height *= self.scale_factor
        width *= self.scale_factor
        self.slit[int(x0 - height//2):int(x0 + height//2), int(y0 - width//2):int(y0 + width//2)] = 1+0j 
        

    def grid(self, width, height, x0, y0, N, offset):
        for n in range(N):
            exp.rectangle(width, height, self.N//2, self.N//2 - n*offset)
            exp.rectangle(width, height, self.N//2, self.N//2 + n*offset)
            exp.rectangle(height, width, self.N//2 - n*offset, self.N//2 )
            exp.rectangle(height, width, self.N//2 + n*offset, self.N//2 )

    def double_slit(self, offset, height, width):
        """ 
        Replicates a double slit experiment
        """
        self.rectangle(width, height, self.size//2, self.size//2 + offset)
        self.rectangle(width, height, self.size//2, self.size//2 - offset)
    def fresnel_rings(self, Z, maxrings):
        """Creates a filter given by the fresnel rings. 
        Params: 
        Z: float
            Length at which the focus is to be found
        """
        radii = [] #A list of the radii
        n = 1 #r0 = 0
        while True:
            rm = n * self.wavelength * np.sqrt( Z/n/self.wavelength + n**2 / 4) #Calculated by pythagoras
            if rm > self.size//2 or n>maxrings:
                break #I want every ring to be completely inside my square
            radii.append(rm)
            n += 1

        for m, rm in enumerate(radii):
            self.slit += (np.sqrt((self.X - self.size//2)**2 + (self.Y - self.size//2)**2) < rm ) * (-1)**(m) * 0.5 + 0.5

    def propagate(self, Z, A):
        """
        Params:
        Z: The distance from the screen to the slit
        A: The Fourier transform A(kx, ky, 0)

        Returns the forier transform A(kx, ky, Z)
        """
        phase = np.exp(1j*Z*np.sqrt( (2*np.pi/self.wavelength)**2 - self.Kx**2 - self.Ky**2))  #Not mine :)
        return A*phase

    def fourier_transform(self):
        """Returns the fourier transform of the slit.
        """
        fourier = fft.fft2(self.slit)
        return fft.fftshift(fourier) #necessary to get nice and tidy graphic results

    def intensity(self, F):
        """Returns the intensity of a field F
        Params:
        F: float
        """
        return np.real(F*np.conjugate(F)) #np real not to get a +0j

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

        intensity = self.intensity(field)
        return intensity/intensity.max()

    def evolve(self, Zmin, Zmax, N=100):
        """Returns a list of the fields for different Z"""
        for Z in np.linspace(Zmin, Zmax, N): #A stupidly complex way of changing Z
            yield (self.field_at_Z(Z), Z)



N = 700 # pixels number
size = 20 * mm # physical length of the experiment (size_x = size_y)
#wavelength = 0.02*mm
wavelength = 700*nm

exp = Monochromatic_Experiment(wavelength, N, size)

"""Double slit experiment"""
a = 1
L = 2000
#exp.rectangle(a, a, size//2, size//2)
#exp.double_slit(1, 40, 0.2)

fig, ax = plt.subplots(2, 1)

ax[1].tick_params(left=False,
                bottom=True,
                labelleft=False,
                labelbottom=True)

#ax.imshow(exp.field_at_Z(1000))
#field = exp.field_at_Z(L)
#field0 = field[N//2, N//2]
Z_fresnel=700
exp.fresnel_rings(Z_fresnel, 100)
#exp.spherical_wave(1000)
ax[1].imshow(np.real(exp.slit))
#ax[0].set_title('Distance: ' +str(int(L)) + ' mm\nSide = ' + str(int(size)) + ' mm')

def animate(Z):
    field = exp.field_at_Z(Z)
    ax[0].cla()
    ax[0].set_title('Distance: ' +str(int(Z)) + ' mm\nSide = ' + str(int(size)) + ' mm, Wavelength = ' + str(wavelength*10**6) + ' nm')
    ax[0].plot(exp.x, field[N//2, :])

    ax[1].cla()
    ax[1].imshow(field)

#ax[0].legend(loc='best')

ani = FuncAnimation(fig, animate, frames=np.linspace(Z_fresnel, Z_fresnel, 1), interval=1, repeat=False)

plt.show()
