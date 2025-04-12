import time


##############################################################################################################################
# Fonctions

def _function_parameter(y, k=None, a = None):
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
        k1 = F(_function_parameter(y))
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
        k1 = F(_function_parameter(y))
        k2 = F(_function_parameter(y, k1, dt))
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
        k1 = F(_function_parameter(y))
        k2 = F(_function_parameter(y, k1, dt/2))
        k3 = F(_function_parameter(y, k2, dt/2))
        k4 = F(_function_parameter(y, k3, dt))
        for i in range(len(y)):
            y[i].append(y[i][-1] + (dt/6)*(k1[i] + 2*k2[i]+ 2*k3[i] + 2*k4[i]))
        t.append(t[-1] + dt)
    print(f"Calculs Termines. Temps de calcul = {round(time.time()-time_init, 4)} s")
    return y, t