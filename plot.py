#import#
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# initalize #
grid = np.array([128,64,32,16])#, 128, 256])
error = []
dsize = []
dsmall = []
mode = "data"

# loop #
for q in range (0, len(grid)):

  #load#
  try:
    data = h5py.File("./"+str(mode)+"/results-"+str(grid[q])+".hdf5", "r")
  except:
    continue 
  
  #load#
  dx_small = data['dx_small'][()][0].T
  dx_big = data['dx_big'][()][0].T
  mass = data['mass'][()][0].T
  r_surf = data['r_surf'][()][0].T
  x2 = data['x2'][()].T
  y2 = data['y2'][()].T
  z2 = data['z2'][()].T
  phi = data['phi'][()].T

  # meshgrid #
  x, y, z = np.meshgrid(x2, y2, z2, indexing="ij")
  rad = np.sqrt(x**2+y**2+z**2)
  v_exact = rad.copy()

  # assign exact potential #
  for i in range (0, rad.shape[0]):
    for j in range (0, rad.shape[1]):
      for k in range (0, rad.shape[2]):
        if(rad[i,j,k] > r_surf):
          v_exact[i,j,k] = - mass/rad[i,j,k]
        else:
          v_exact[i,j,k] = - 0.5*mass*(3*r_surf**2 - rad[i,j,k]**2)/(r_surf**3)

 
  # set vmin #
  cmap = "jet"
  vmin = v_exact.min()
  vmax = v_exact.max()

  #plot#
  if(q == 0):
    fig, axs = plt.subplots()
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    tick = np.linspace(vmin, vmax, 10)
    plt.contourf(x2, y2, phi[:,:,phi.shape[2]//2], 100, norm=norm, cmap=cmap)
    plt.gca().set_aspect('equal')
    plt.ylabel(r'$y-direction$',size=15)
    plt.xlabel(r'$x-direction$',size=15)
    cbar_ax = fig.add_axes([0.86, 0.15, 0.03, 0.7])
    cbar = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, ticks=tick)
    plt.savefig('numerical.png')
    plt.clf()
    plt.close()

    fig, axs = plt.subplots()
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    tick = np.linspace(vmin, vmax, 10)
    plt.contourf(x2, y2, v_exact[:,:,phi.shape[2]//2], 100, norm=norm, cmap=cmap)
    plt.gca().set_aspect('equal')
    plt.ylabel(r'$y-direction$',size=15)
    plt.xlabel(r'$x-direction$',size=15)
    cbar_ax = fig.add_axes([0.86, 0.15, 0.03, 0.7])
    cbar = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, ticks=tick)
    plt.savefig('exact.png')
    plt.clf()
    plt.close()

    vdiff = np.log10(np.abs(phi-v_exact))
    vmin = vdiff.min()
    vmax = vdiff.max()
    fig, axs = plt.subplots()
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    tick = np.linspace(vmin, vmax, 10)
    plt.contourf(x2, y2, vdiff[:,:,vdiff.shape[2]//2], 100, norm=norm, cmap=cmap)
    plt.gca().set_aspect('equal')
    plt.ylabel(r'$y-direction$',size=15)
    plt.xlabel(r'$x-direction$',size=15)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    cbar = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, ticks=tick)
    plt.savefig('absdiff.png')
    plt.clf()
    plt.close()

  # l2-norm #
  norm = v_exact - phi
  norm = norm**2
  norm = norm.flatten()
  norm = np.sum(norm)/len(norm)
  error.append(norm)
  dsize.append(dx_big)
  dsmall.append(dx_small)

error = np.array(error)
dsize = np.array(dsize)
dsmall = np.array(dsmall)
plt.title("MSE Numerical vs Exact")
plt.scatter(dsize, error)
plt.plot(dsize, error, linestyle=":", label=r"$\Delta h$")
plt.scatter(dsmall, error)
plt.plot(dsmall, error, linestyle=":", label=r"$\Delta hr^{N/2-1}$")
plt.grid()
plt.xlabel(r"Step Size")
plt.ylabel(r"Mean Sqaure Error")
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.savefig("convergence"+str(mode)+".png")
