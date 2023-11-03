import matplotlib.pyplot as plt
import numpy as np
import time
import random


Niter = int(1e4)
precision = 1
dt = 10**(-precision)
Nbloc = 50

#distance entre les blocs
l0 = 1/(Nbloc)

kc = 1
vc = 0.001
kp = 10
m = 1

mu_s = 0.25
mu_d = 0.2
dc = 1e-4


T = np.arange(0, (Niter+1)*dt, dt)
U = np.zeros((Niter+1,Nbloc))
V = np.zeros((Niter+1,Nbloc))
A = np.zeros((Niter+1,Nbloc))

liste_Gliss = np.zeros((Niter+1,Nbloc))


for i in range(Nbloc):
    U[0,i] = (i+(random.randint(100-20,100+20)/100))/Nbloc

T0 = [0 for i in range(Nbloc)]


def frottement(t,p):
    global U, T0
    if (U[t,p] - U[T0[p], p] < dc):
        return (U[t,p]-U[T0[p],p])*(mu_d-mu_s)/dc + mu_s
    else:
        return mu_d
    
def f_ressort(t,p):
    global T, U
    if p == 0:
        return kc*(p*1/(Nbloc) + vc*T[t] - U[t,p]) + kp*(U[t,p+1] - U[t,p] - l0)
    
    elif p == Nbloc-1:
        return kc*(p*1/(Nbloc) + vc*T[t] - U[t,p]) - kp*(U[t,p] - U[t,p-1] - l0) 
    
    else :
        return kc*(p*1/(Nbloc) + vc*T[t] - U[t,p]) + kp*(U[t,p+1] - 2*U[t,p] + U[t,p-1])


start = time.time()
for t in range(Niter):
    for p in range(Nbloc):  
        if liste_Gliss[t, p]:
            U[t+1,p] = U[t,p] + dt*V[t,p] + dt**2*(f_ressort(t,p)  - frottement(t,p))/(2*m)
            # U[t+1,p] = 2*U[t,p] - U[t-1,p] + dt**2*(f_ressort(t,p)  - frottement(t,p))/m
            # V[t+1,p] = (U[t+1,p]-U[t,p])/dt
            
        else:
            U[t+1,p] = U[t,p]
       
    for p in range(Nbloc):
        if liste_Gliss[t, p]:
            V[t+1,p] = V[t,p] + dt*(f_ressort(t, p)  - frottement(t,p) + f_ressort(t+1,p)- frottement(t+1,p))/(2*m)
        
        if f_ressort(t+1, p) >= mu_s:
            liste_Gliss[t+1, p] = 1
            T0[p] = t+1
        

        if V[t+1,p] < 0 :
            liste_Gliss[t+1, p] = 0

        elif liste_Gliss[t, p]:
            liste_Gliss[t+1, p] = 1



end = time.time()

print(f"temps d'execution : {end-start:.3}")     


# plt.plot(T,U[:,0],label=f"bloc {p}")
for p in range(Nbloc):
    plt.plot(T,U[:,p],label=f"bloc {p}")
    # plt.legend()
    plt.title(f"mu_d = {mu_d}, mu_s = {mu_s}, dc = {dc} ")
    plt.xlabel('temps t')
    plt.ylabel('déplacement u')

    # plt.xlim(140,160)
    # plt.ylim(0.4,0.6)

plt.show()


"""
Regarder bloc par bloc pour aléatoire
mettre force ressort à l'equilibre au début'

"""

start = time.time()

with open('test.json', 'w', newline='') as fichier:

    fichier.write('{\n')
    for i in range(Niter):
        ligne = f"\"{T[i]:.{precision}F}\" : {[list(U[i,:]), list(liste_Gliss[i,:])]}" + ', \n'
        fichier.write(ligne)

    der_ligne = f"\"{T[Niter]:.{precision}F}\" : {[list(U[Niter,:]), list(liste_Gliss[Niter,:])]}" + '\n'
    fichier.write(der_ligne)

    fichier.write('}')

end = time.time()
print(f"temps d'écriture : {end-start:.3}") 

# with open('test.json', 'r', newline='') as jsonfile:
#     data = json.load(jsonfile)

#     print(data["0.00"])
#     print(data["500.00"])

