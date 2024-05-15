from vpython import *
import matplotlib.pyplot as plt

# costanti
G = 6.67e-11
Ms = 1e30
AU = 1.49e11
Rs = 190e7

# masse
mA = Ms
mB = 0.8 * mA
mC = 0.4 * mA

# raggi
rA = 1.2 * Rs
rB = rA * pow(mB / mA, 3)
rC = 5 * rA * pow(mC / mA, 3)

# distanze dal CM
xA = 0.1 * AU
zC = 0.6 * AU
xB = xA * mA / mB
zA = zB = zC * mC / (mA + mB)

# spheres
sA = sphere(pos=vector(xA, 0, -zA), radius=rA, color=color.yellow, make_trail=True)
sB = sphere(pos=vector(-xB, 0, -zB), radius=rB, color=color.cyan, make_trail=True)
sC = sphere(pos=vector(0, 0, zC), radius=rC, color=color.magenta, make_trail=True)
vA = pow(G * mB * xA, 0.5) / (xA + xB)

# momenti t=0
sA.p = mA * vA * vector(0, 1, 0)
sB.p = -sA.p
sC.p = vector(0, 0, 0)

t = 0
dt = 1000
r = 400
t_sec = 30

# Style info
pot_label = 'potenziale gravitazionale'
kin_label = 'cinetica'
mec_label = 'meccanica'
x_label = 'Tempo (s)'
y_label = 'Energia (J)'

# VPython plot
gr = graph(xtitle=x_label, ytitle=y_label, ymin=-16e39, ymax=18e39)
pot_en_graph = gcurve(graph=gr, color=color.blue, label=f'Energia {pot_label}')
kin_en_graph = gcurve(graph=gr, color=color.red, label=f'Energia {kin_label}')
mec_en_graph = gcurve(graph=gr, color=color.black, label=f'Energia {mec_label}')

# matplotlib plot
pot_en_data = []
kin_en_data = []
mec_en_data = []
time_frames = [ti * dt for ti in range(r * t_sec)]


def calc_pot_energy(rAB, rBC, rCA):
    pot_energy = -G * (mA * mB / mag(rAB) + mB * mC / mag(rBC) + mC * mA / mag(rCA))
    return pot_energy


def calc_kin_energy(pA, pB, pC):
    kin_energy = 0.5 * (mA * pow(pA / mA, 2) + mB * pow(pB / mB, 2) + mC * pow(pC / mC, 2))
    return kin_energy


def plot_data(time_val, data, lab, col, ax, legend=False):
    ax.plot(time_val, data, linestyle='-', color=col, label=f'Energia {lab}')
    ax.grid(True)
    if legend: ax.legend()
    # ax.set_title(lab)
    # ax.xlabel(x_label)
    # ax.ylabel(y_label)
    # plt.title('')


sleep(1)
while t < t_sec * r * dt:
    rate(r)
    rAB = sB.pos - sA.pos
    rBC = sC.pos - sB.pos
    rCA = sA.pos - sC.pos
    FBA = G * mA * mB * norm(rAB) / mag(rAB) ** 2
    FCB = G * mC * mB * norm(rBC) / mag(rBC) ** 2
    FAC = G * mA * mC * norm(rCA) / mag(rCA) ** 2
    sA.F = FBA - FAC
    sB.F = -FBA + FCB
    sC.F = -FCB + FAC
    sA.p = sA.p + sA.F * dt
    sB.p = sB.p + sB.F * dt
    sC.p = sC.p + sC.F * dt
    sA.pos = sA.pos + sA.p * dt / mA
    sB.pos = sB.pos + sB.p * dt / mB
    sC.pos = sC.pos + sC.p * dt / mC
    t = t + dt

    pot_en = calc_pot_energy(rAB, rBC, rCA)
    kin_en = calc_kin_energy(mag(sA.p), mag(sB.p), mag(sC.p))
    mec_en = pot_en + kin_en
    pot_en_graph.plot(t, pot_en)
    kin_en_graph.plot(t, kin_en)
    mec_en_graph.plot(t, mec_en)

    pot_en_data.append(pot_en)
    kin_en_data.append(kin_en)
    mec_en_data.append(mec_en)

fig, (ax1, ax2) = plt.subplots(2)
plot_data(time_frames, kin_en_data, kin_label, 'r', ax1)
plot_data(time_frames, mec_en_data, mec_label, 'black', ax1)
plot_data(time_frames, pot_en_data, pot_label, 'b', ax1)

plot_data(time_frames, mec_en_data, mec_label, 'black', ax2)
input('Mostra grafico...')

for ax in fig.get_axes():
    ax.label_outer()
plt.show()


