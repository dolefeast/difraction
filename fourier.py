import matplotlib.pyplot as plt
import scipy.fft as fft
import numpy as np

N = 200
height = 80 
width = 80
fig, ax = plt.subplots(1, 2)

m = 1
cm = m/10
mm = cm/10
nm = 1e-9

l_ambda = 0.0001
Z = 40
h=1
A=5

xfig = np.zeros((N, N))


xfig[(N-height)//2:(N+height)//2, (N-width)//2:(N+width)//2] += 1

ax[0].imshow(xfig)

X = np.arange(0, N, 1)
Y = np.arange(0, N, 1)
U, V = np.meshgrid(X, Y)

yfig = fft.fft2(xfig)
yfigshift = fft.fftshift(yfig)

for Z in range(200):
    if Z == 40:
        print('yiaier')
    phase = np.exp(2*np.pi*1j*(Z/h)*np.sqrt((1/l_ambda)**2-(U/A)**2-(V/A)**2))#print(phase)

    field = yfigshift*phase
    #print(field[0,0], yfigshift[0,0])

    U = fft.ifft2(yfigshift*phase)

    ax[1].cla()
    ax[1].imshow(np.real(np.conjugate(U)*U))
    plt.pause(0.001)
plt.savefig('img.pdf')
plt.show()


