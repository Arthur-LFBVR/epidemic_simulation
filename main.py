# Importation des modules
import numpy as np
import matplotlib.pyplot as plt
import time

##############################################################################################################################
# Paramètres initiaux

S_init = 999    # Population initialement saine
I_init = 1      # Population initialement malade
R_init = 0      # Population initialement gérie
D_init = 0

pop_tot = S_init + I_init + R_init  # Population totale

p = 0.25        # Paramètre de tranmission
alpha = 0.015    # Paramètre de guérison
mu = 0          # Paramètre de mortalité

T = 300             # Durée de la simulation
dt = 0.001          # Intervale de temps

##############################################################################################################################
# Fonctions

def model_SIR(y : list):
    dR = alpha * y[1]
    dS = -p * y[1] * y[0]
    dI = -dS - dR
    return [dS, dI, dR]


def model_SIRD(y : list):
    dR = alpha * y[1]
    dD = mu * y[1]
    dS = -p * y[1] * y[0]
    dI = -dS - dR - dD
    return [dS, dI, dR, dD]


def function_parameter(y, k=None, a = None):
    temp = []
    if k == None and a == None:
        for i in range(len(y)):
            temp.append(y[i][-1])
    else :
        for i in range(len(y)):
            temp.append(y[i][-1] + a*k[i])
    return temp


def euler(F, y : list, T : float, dt : float):
    t = [0]     # Initialisation du temps
    time_init = time.time()
    while t[-1] <= T :
        print(f"Avancement des calculs : {round((t[-1]/T)*100, 4)} %")
        k1 = F(function_parameter(y))
        for i in range(len(y)):
            y[i].append(y[i][-1] + k1[i] * dt)
        t.append(t[-1] + dt)
    print(f"Calculs Termines. Temps de calcul = {round(time.time()-time_init, 4)} s")
    return y, t


def RK2(F, y : list, T : float, dt : float):
    t = [0]     # Initialisation du temps
    time_init = time.time()
    while t[-1] <= T :
        print(f"Avancement des calculs : {round((t[-1]/T)*100, 4)} %")
        k1 = F(function_parameter(y))
        k2 = F(function_parameter(y, k1, dt))
        for i in range(len(y)):
            y[i].append(y[i][-1] + (dt/2)*(k1[i] + k2[i]))
        t.append(t[-1] + dt)
    print(f"Calculs Termines. Temps de calcul = {round(time.time()-time_init, 4)} s")
    return y, t


def RK4(F, y : list, T : float, dt : float):
    t = [0]     # Initialisation du temps
    time_init = time.time()
    while t[-1] <= T :
        print(f"Avancement des calculs : {round((t[-1]/T)*100, 4)} %")
        k1 = F(function_parameter(y))
        k2 = F(function_parameter(y, k1, dt/2))
        k3 = F(function_parameter(y, k2, dt/2))
        k4 = F(function_parameter(y, k3, dt))
        for i in range(len(y)):
            y[i].append(y[i][-1] + (dt/6)*(k1[i] + 2*k2[i]+ 2*k3[i] + 2*k4[i]))
        t.append(t[-1] + dt)
    print(f"Calculs Termines. Temps de calcul = {round(time.time()-time_init, 4)} s")
    return y, t


##############################################################################################################################
# Programme principal

y, t = euler(model_SIRD, [[S_init/pop_tot], [I_init/pop_tot], [R_init/pop_tot]], T, dt)

fig, ax = plt.subplots()
fig.canvas.manager.set_window_title("Epidemic Simulation")

plt.xlabel("Temps")
plt.ylabel("Proportion de la population")
plt.title("Évolution de l’épidémie au cours du temps")

ax.plot(t, y[0], label = "Population saine", color = "blue")
ax.plot(t, y[1], label = "Population malade", color = "green")
ax.plot(t, y[2], label = "Population guérie", color = "red")

plt.tight_layout()
plt.show()