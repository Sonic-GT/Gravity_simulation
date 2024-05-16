from vpython import *
import matplotlib.pyplot as plt

# costanti
G = 6.67e-11
Ms = 1e30
AU = 1.49e11
Rs = 190e7

# masse
mA = Ms
mB = 0.5 * mA
mBody = 1 / (1 / mA + 1 / mB)
mCM = mA + mB

# raggi
rBody = 1.2 * Rs
rCM = 0.3 * rBody

# distanze dal CM
xA = 0.08 * AU
xB = xA * mA / mB
xBody = xA + xB

# spheres
sBody = sphere(pos=vector(xBody, 0, 0), radius=rBody, color=color.yellow, make_trail=True)
sCM = sphere(pos=vector(0, 0, 0), radius=rCM, color=color.cyan, make_trail=True)
vA = 0.4 * pow(G * mB * xA, 0.5) / (xA + xB)
vB = mA / mB * vA

# quantità di moto t=0
sBody.p = mBody * (vA + vB) * vector(0,1,0)
sCM.p = 0

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


def calc_pot_energy(rBody_CM):
    pot_energy = -G * mBody * mCM / mag(rBody_CM)
    return pot_energy


def calc_kin_energy(pBody):
    kin_energy = 0.5 * mBody * pow(pBody / mBody, 2)
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
    rBody_CM = sCM.pos - sBody.pos
    FCM_Body = G * mA * mB * norm(rBody_CM) / mag(rBody_CM) ** 2
    sBody.F = FCM_Body
    sBody.p = sBody.p + sBody.F * dt
    sBody.pos = sBody.pos + sBody.p * dt / mBody
    t = t + dt

    pot_en = calc_pot_energy(rBody_CM)
    kin_en = calc_kin_energy(mag(sBody.p))
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
