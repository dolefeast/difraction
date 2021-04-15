from fourier import *

size = 20 * mm
wavelength = 900 * nm
N = 2000

exp = Monochromatic_Experiment(wavelength, N, size)

a = 0.8
exp.double_slit(a, 20, 0.4)

fig, (ax1, ax2) = plt.subplots(2, 1)


Z = 2000
#F = a**2/Z/wavelength
ax1.set_title("Distancia={} mm, L. onda={} nm\nTamaño={} mm ".format(Z, wavelength * 10**6 , size))

field = exp.field_at_Z(Z)

ax1.plot(exp.x, field[N//2, :])
ax2.imshow(field)

plt.show()
fig.savefig('./fig/double_slit.pdf')
