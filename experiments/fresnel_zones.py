from fourier import *

wavelength = 700 * nm
N = 700
size = 20 * mm
exp = Monochromatic_Experiment(wavelength, N, size)

Z = 700
exp.fresnel_rings(Z, 100, sign=0)

fig, (ax1, ax2) = plt.subplots(2, 1)

field = exp.field_at_Z(Z)
ax1.plot(exp.x, field[:, N//2])
ax2.imshow(field)

plt.savefig('./fig/fresnel_rings.pdf')
plt.show()
