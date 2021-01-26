#!/usr/bin/env pvpython
import numpy as np
from scipy import linalg
from vtk.numpy_interface import dataset_adapter as da

s_xx = inputs[0].PointData['stress_xx']
s_yy = inputs[0].PointData['stress_yy']
s_zz = inputs[0].PointData['stress_zz']
s_xy = inputs[0].PointData['stress_xy']
s_yz = inputs[0].PointData['stress_yz']
s_xz = inputs[0].PointData['stress_xz']

direction = []
sigma_d = []

for i in range(len(s_xx)):
   stress = np.array ( [[s_xx[i], s_xy[i], \
  s_xz[i]], [s_xy[i], s_yy[i], s_yz[i] ], [s_xz[i], \
  s_yz[i], s_zz[i] ]] )
   u,v  = np.linalg.eig(stress)
   indx = np.argsort(u)
   sigma_d.append(u[indx[2]] - u[indx[0]]) 
   direction.append(v[indx[2]].tolist())

vtk_arr = da.VTKArray(direction)
output.PointData.append(vtk_arr, "direction_1")


####### routine for stress rotation and coloumb stress calculation


for j in range(0, 90, 10):

  theta = np.radians(j)
  coloumb = []
  phi = np.radians(40) ; # along z

  # rotation matrix for rotation about z
  R = [[np.cos(theta), -np.sin(theta), 0 ], [np.sin(theta), np.cos(theta), 0 ], [0, 0, 1]]

  Rx = [[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi) ]]

  for i in range(len(s_xx)):
    stress = np.array ( [[s_xx[i], s_xy[i], s_xz[i]], [ s_xy[i], \
      s_yy[i], s_yz[i] ], [s_xz[i], s_yz[i], s_zz[i] ] ] )
    
    stress_dash_x = np.dot (np.dot (np.transpose(Rx), stress), Rx)

    stress_dash = np.dot (np.dot (np.transpose(R), stress_dash_x), R) # fault plane
    sigma_nn = stress_dash[1][1] # normal stress sigma_yy 
    tau = stress_dash[0][1] # assuming strike-slip component sigma_xy
    mu  = 0.6 # coefficient of friction
    cff = tau + mu*sigma_nn
    coloumb.append(cff)
  output.PointData.append(coloumb, "coloumb_%s" %j)

output.PointData.append(sigma_d, "sigma_d")
