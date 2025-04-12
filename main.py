# Importation des modules
import matplotlib.pyplot as plt
import solv_diff_eq

##############################################################################################################################
# Paramètres initiaux

S_init = 999    # Population initialement susceptible
I_init = 1      # Population initialement infectée
R_init = 0      # Population initialement rétablie
D_init = 0      # Population initialement décédée

pop_tot = S_init + I_init + R_init  # Population totale

p = 0.25            # Paramètre de tranmission
alpha = 0.015       # Paramètre de guérison
mu = 0.001          # Paramètre de mortalité

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

def model_SIS(y : list):
    dS = -p * y[0] * y[1] + alpha * y[1]
    dI = -dS
    return [dS, dI]



##############################################################################################################################
# Programme principal

info_graph = [{"label" : "Population susceptible", "color" : "blue"},
              {"label" : "Population infectée", "color" : "green"},
              {"label" : "Population rétablie", "color" : "red"},
              {"label" : "Population décédée", "color" : "black"},]

y, t = solv_diff_eq.RK4(model_SIRD, [[S_init/pop_tot], [I_init/pop_tot], [R_init/pop_tot], [D_init/pop_tot]], T, dt)

fig, ax = plt.subplots(figsize=(10, 5))
fig.canvas.manager.set_window_title("Epidemic Simulation")

plt.xlabel("Temps")
plt.ylabel("Proportion de population")
plt.title("Évolution de l’épidémie au cours du temps")

for i in range(len(y)):
    ax.plot(t, y[i], label = info_graph[i]["label"], color = info_graph[i]["color"])

ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)

plt.tight_layout()
plt.show()