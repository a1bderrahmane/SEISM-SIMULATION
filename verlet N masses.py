import matplotlib.pyplot as plt
import numpy as np
import time
import random
import json


Niter = int(1e5)
precision = 1
dt = 10**(-precision)
Nbloc = 500

alea = 20               #pourcentage d'aleatoire sur la position initiale

L = 1                   #largeur de la faille
# a = l0 = L/Nbloc      #distance entre les blocs a l'equilibre
a = l0 = 0.01

kc = 1
vc = 0.0001
kp = 100
m = 1

mu_s = 0.25
mu_d = 0.2
dc = 1e-6




# T = [i*dt for i in range(Niter+1)]
T = np.arange(0, (Niter+1)*dt, dt) 

U = np.zeros((Niter+1,Nbloc))
V = np.zeros((Niter+1,Nbloc))
A = np.zeros((Niter+1,Nbloc))

# liste_temp = [False for i in range(Nbloc)]
liste_Gliss = np.zeros((Niter+1,Nbloc), dtype = bool)

U[0,:] = [i*a + a*(random.randint(100-alea,100+alea)/100) for i in range(Nbloc)]
# U[1,:] = U[0,:]
print(U[0,:])

T0 = [0 for i in range(Nbloc)]


def frottement(t,p):
    global U, T0
    if (U[t,p] - U[T0[p],p]) < dc:
        return (U[t,p]-U[T0[p],p])*(mu_d-mu_s)/dc + mu_s
    else:
        return mu_d
    
def f_ressorts(t,p):
    global T, U
    if p == 0:
        return kc*(p*a + vc*T[t] - U[t,p]) + kp*(U[t,p+1] - U[t,p] - l0)
    
    elif p == Nbloc-1:
        return kc*(p*a + vc*T[t] - U[t,p]) - kp*(U[t,p] - U[t,p-1] - l0) 
    
    else :
        return kc*(p*a + vc*T[t] - U[t,p]) + kp*(U[t,p+1] - 2*U[t,p] + U[t,p-1])


start = time.time()
for t in range(Niter):
    for p in range(Nbloc):  
        if liste_Gliss[t,p]:
            U[t+1,p] = U[t,p] + dt*V[t,p] + dt**2*(f_ressorts(t,p)  - frottement(t,p))/(2*m)
            
        else:
            U[t+1,p] = U[t,p]
        
    for p in range(Nbloc):
        if liste_Gliss[t,p]:
            V[t+1,p] = V[t,p] + dt*(f_ressorts(t, p)  - frottement(t,p) + f_ressorts(t+1,p)- frottement(t+1,p))/(2*m)
        
        if f_ressorts(t+1, p) >= mu_s:
            liste_Gliss[t+1,p] = True        #declenche le glissement
            T0[p] = t+1
        
        if V[t+1,p] < 0 :
            liste_Gliss[t+1,p] = False
        elif liste_Gliss[t,p]:              #si on glissait au temps d'avant, et que la vitesse n'est pas négative, on continue
            liste_Gliss[t+1,p] = True
            
    
    # liste_Gliss[t+1,:] = liste_temp
end = time.time()

print(f"temps d'execution : {end-start:.3}")     










start = time.time()
for p in range(Nbloc):
    plt.plot(T,U[:,p],label=f"bloc {p}")
    # plt.xlim(0,500)
    # plt.ylim(0,1.2)
    
# plt.legend()
plt.xlabel('Temps t')
plt.ylabel('Position u')
plt.title(f"$\mu_s$ = {mu_s}, $\mu_d$ = {mu_d}, $d_c$ = {dc}\n $v_c$ = {vc}, $k_c$ = {kc}, $k_p$ = {kp}")
plt.show()
end = time.time()
print(f"temps d'affichage : {end-start:.3}") 









dicoEvenement = {}
i = 0
while i < Niter:
    if True in liste_Gliss[i,:] :
        debut_event = i
        while True in liste_Gliss[i,:] and i < Niter :
            i += 1
        fin_event = i
        
        num_event = len(dicoEvenement)
        dicoEvenement[f"event_{num_event}"] = [debut_event, fin_event, U[debut_event,:], U[fin_event,:]]
    i += 1
    
plt.clf()
for (key, value) in dicoEvenement.items():
    pos_debut_event = value[2]
    pos_fin_event = value[3]
    deplacement = pos_fin_event - U[0,:]
    liste_x = [i for i in range(Nbloc)]
    plt.plot(liste_x, deplacement, color='k', lw=0.5)
plt.show()






"""
Calcul magnitude
"""
liste_magnitude = []
for event in dicoEvenement.values():
    t_debut = event[0]
    t_fin = event[1]
    pos_debut = event[2]
    pos_fin = event[3]


    aire_deb = 0
    aire_fin = 0
    for num_bloc in range(len(pos_debut)-1):
        pos_bloc_n = (pos_debut[ num_bloc]  + pos_fin[num_bloc])/2
        pos_bloc_n_plus_1 = (pos_debut[num_bloc+1]  + pos_fin[num_bloc+1])/2
        largeur = pos_bloc_n_plus_1 - pos_bloc_n
        
        #Au debut de l'event : 
        hauteur_deb = (pos_debut[num_bloc] + pos_debut[num_bloc+1])/2
        aire_deb += largeur * hauteur_deb
        
        #A la fin de l'event : 
        hauteur_fin = (pos_fin[num_bloc] + pos_fin[num_bloc+1])/2
        aire_fin += largeur * hauteur_fin

    magnitude = aire_fin - aire_deb
        
    liste_magnitude.append(magnitude)

liste_magnitude.sort()

liste_freq = []
N_event = len(liste_magnitude)
for i in range(1,N_event):
    liste_freq.append(1 - i/N_event)

plt.clf()
plt.plot(liste_magnitude, liste_freq, ' x')
plt.xlabel("magnitude")
plt.ylabel("frequence")
# plt.semilogy()
plt.loglog()
plt.show()
    





"""
start = time.time()

with open('test.json', 'w', newline='') as fichier:

    fichier.write('{\n')

    for i in range(Niter):
        ligne = f"\"{T[i]:.{precision}F}\" : {[list(U[i,:]), list(liste_Gliss[i,:].astype(float))]}" + ', \n'
        fichier.write(ligne)

    der_ligne = f"\"{T[Niter]:.{precision}F}\" : {[list(U[Niter,:]), list(liste_Gliss[Niter,:].astype(float))]}" + '\n'
    fichier.write(der_ligne)

    fichier.write('}')

end = time.time()
print(f"temps d'écriture : {end-start:.3}") 

# with open('test.json', 'r', newline='') as jsonfile:
#     data = json.load(jsonfile)

#     print(data["0.00"])
#     print(data["500.00"])


"""




















