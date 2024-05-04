from vpython import *
import matplotlib.pyplot as plt

centre_mB = True  # Centra la scena in mB?

# costanti
G = 6.67e-11
Ms = 1e30
AU = 1.49e11
Rs = 190e7

# masse
mA = Ms
mB = 0.5 * mA

# raggi
rA = 1.2 * Rs
rB = 3 * rA * pow(mB / mA, 3)

# distanze dal CM
xA = 0.08 * AU
xB = xA * mA / mB

# spheres
sA = sphere(pos=vector(xA, 0, 0), radius=rA, color=color.yellow, make_trail=not centre_mB)
sB = sphere(pos=vector(-xB, 0, 0), radius=rB, color=color.cyan, make_trail=not centre_mB)
vA = 0.4 * pow(G * mB * xA, 0.5) / (xA + xB)

if centre_mB: sB_trail_curve = curve(color=color.yellow)

# momenti t=0
sA.p = mA * vA * vector(0, 1, 0)
sB.p = -sA.p

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


def calc_pot_energy(rAB):
    pot_energy = -G * mA * mB / mag(rAB)
    return pot_energy


def calc_kin_energy(pA, pB):
    kin_energy = 0.5 * (mA * pow(pA / mA, 2) + mB * pow(pB / mB, 2) )
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
    FBA = G * mA * mB * norm(rAB) / mag(rAB) ** 2
    sA.F = FBA
    sB.F = -FBA
    sA.p = sA.p + sA.F * dt
    sB.p = sB.p + sB.F * dt
    sA.pos = sA.pos + sA.p * dt / mA
    sB.pos = sB.pos + sB.p * dt / mB
    t = t + dt

    if centre_mB:
        scene.center = sB.pos
        sB_trail_curve.origin = scene.center
        pos_from_camera = sA.pos - sB.pos
        sB_trail_curve.append(pos=pos_from_camera)

    pot_en = calc_pot_energy(rAB)
    kin_en = calc_kin_energy(mag(sA.p), mag(sB.p))
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
