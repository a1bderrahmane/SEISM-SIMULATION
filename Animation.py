import json
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.pyplot as plt
 
plt.rcParams['animation.ffmpeg_path'] = r'C:/Users/mmosca1/Downloads/ffmpeg-master-latest-win64-gpl-shared/bin/ffmpeg.exe'



with open('test.json', 'r', newline='') as jsonfile:
    data = json.load(jsonfile)

    # print(data["0.00"])
    # print(data["500.00"])


Nbloc = len(data[list(data.keys())[0]][0])
Niter = len(data.keys()) - 1
dt = float(list(data.keys())[1]) - float(list(data.keys())[0])


# print(Nbloc)
# print(Niter)
# print(dt)


U = np.zeros((Niter+1,Nbloc))
T = np.arange(0, (Niter+1)*dt, dt)

j = 0
for cle in data.keys():
    U[j,:] = data[cle][0]
    j += 1

# print(type(U[0]))
# print(T[0])
# print(list(U[0]))

reduction = 10
U_reduit = np.zeros((int((Niter+1)/reduction),Nbloc))
T_reduit = np.arange(0, int((Niter+1)*dt/reduction), dt)

# print(len(U_reduit))
for i in range(0,len(U),reduction):
    U_reduit[int(i/reduction)-1,:] = U[i,:]



# print(U_reduit[-1])
# print(U[-1])



fig, ax = plt.subplots(figsize=(19.20, 10.80))
# échelle de temps
ax.set_xlim(0, int((Niter+1)*dt/reduction)-1)
# échelle de distance
ax.set_ylim(0, U[-1,-1]+0.001)

#  Initialiser les blocs
blocs = [ax.plot([], [],'o', lw=20)[0] for i in range(Nbloc)]
# bloc, = ax.plot([], [], 'o', lw=20)

# Initialiser les courbes
lines = [ax.plot([], [], label=f'Courbe {i+1}')[0] for i in range(Nbloc)]



def update_courbes(i): 
    # Mettre à jour les lignes avec les nouvelles données
    for bloc in range(Nbloc):
        x_partial_1 = T_reduit[:i+1]
        y_partial_1 = U_reduit[:i+1, bloc]
        lines[bloc].set_data(x_partial_1, y_partial_1)

        x_partial_2 = [T_reduit[i]]
        y_partial_2 = [U_reduit[i, bloc]]
        blocs[bloc].set_data(x_partial_2, y_partial_2)

    return lines, blocs

anim = animation.FuncAnimation(fig, update_courbes, interval=1, frames=int(Niter/reduction))

start = time.time()

writer = animation.FFMpegWriter(fps=24, bitrate=5000)
anim.save("test.mp4", writer = writer, dpi=150)

end = time.time()
print(f"temps d'enregistrement : {end-start:.3}")   

# plt.show()

print("ok")