from fourier import *

size = 100 * mm
wavelength = 9000 * nm
N = 2000

exp = Monochromatic_Experiment(wavelength, N, size)

a = 1
exp.rectangle(a, a, size//2, size//2)

fig, (ax1, ax2) = plt.subplots(2, 1)


Z = 1000
F = a**2/Z/wavelength
ax1.set_title("Distancia={} mm, L. onda={} mm, Tamaño={} mm\nF = {} ".format(Z, wavelength , size, round(F, 3)))

field = exp.field_at_Z(Z)

ax1.plot(exp.x, field[:, N//2])
ax2.imshow(field)

plt.show()
fig.savefig('./fraunhofer.pdf')
