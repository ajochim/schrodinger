#!/usr/bin/env python3
""""""
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.interpolate import interp1d
from scipy import integrate
import matplotlib.pyplot as plt

# Read given parameters of "schrodinger.inp":
schrodingerInp = open("schrodinger.inp", "r")
schrodingerLines = schrodingerInp.readlines()
# print(schrodingerLines)

massEntry = schrodingerLines[0].split()[0]
mass = float(massEntry)
# print("\nmass m:\n", mass)

xMinEntry = schrodingerLines[1].split()[0]
xMin = float(xMinEntry)
# print("\nxMin:\n", xMin)

xMaxEntry = schrodingerLines[1].split()[1]
xMax = float(xMaxEntry)
# print("\nxMax:\n", xMax)

nPointEntry = schrodingerLines[1].split()[2]
nPoint = int(nPointEntry)
# print("\nPoint:\n", nPoint)

firstEigvEntry = schrodingerLines[2].split()[0]
firstEigv = float(firstEigvEntry)
# print("\nfirst eigenvalue to print:\n", firstEigv)

lastEigvEntry = schrodingerLines[2].split()[1]
lastEigv = float(lastEigvEntry)
# print("\nsecond eigenvalue to print:\n", secondEigv)

interType = schrodingerLines[3].split()[0]
# print("\ninterpolation type:\n", interType)

nrInterPointsEntry = schrodingerLines[4].split()[0]
nrInterPoints = int(nrInterPointsEntry)
# print("\nnr. of interpolation points:\n", nrInterPoints)

potLen = len(schrodingerLines) - 5
# print("\nnr. of pot-values:\n", potLen)
xPot = np.zeros(potLen)
yPot = np.zeros(potLen)
for ii in range(5, potLen + 5):
    xPot[ii - 5] = float(schrodingerLines[ii].split()[0])
    yPot[ii - 5] = float(schrodingerLines[ii].split()[1])
print("\nx pot-values:\n", xPot, "\n\ny pot-values:\n", yPot)
schrodingerInp.close()

# Interpolation of the potential:
if interType == "linear":
    print("\nlinear interpolation selected")
    xinterp = np.linspace(xMin, xMax, nPoint)
    print("\ngenerated x values\n", xinterp)
    yinterp = np.interp(xinterp, xPot, yPot)
    print("\ngenerated y values\n", yinterp)

elif interType == "cspline":
    print("\ncspline interpolation selected")
    fCspline = interp1d(xPot, yPot, kind='cubic')
    xinterp = np.linspace(xMin, xMax, nPoint)
    yinterp = fCspline(xinterp)

elif interType == "polynomial":
    polDegree = 2
    print("\npolynomial interpolation of degree", polDegree, "selected")
    coeffPoly = np.polyfit(xPot, yPot, polDegree)  # minimizes the squared error
    fPoly = np.poly1d(coeffPoly)  # read coeffizients
    xinterp = np.linspace(xMin, xMax, nPoint)
    yinterp = fPoly(xinterp)
else:
    print("error: invalid interpolation method")

# write obtained points to file "potential.dat":
InterpPot = np.hstack((xinterp.reshape((-1, 1)), yinterp.reshape((-1, 1))))
np.savetxt("potential.dat", InterpPot)

# create matrix to calculate eigenvalues and eigenvectors:
delta = abs(xinterp[1]-xinterp[0])
print("\nlength of xinterp:\n", len(xinterp), "\n\nobtained delta:\n", delta)
a = 1/(mass*(delta)**2)
ludiag = np.zeros(nPoint - 1)
ludiag[:] = (-a)/2
print("\nlower and upper diagonal:\n", ludiag)
mainDiag = np.zeros(nPoint)

for jj in range(0, nPoint):
    mainDiag[jj] = a + yinterp[jj]

print("\nmainDiag\n", mainDiag)
#diags = np.array([0, -1, 1])
#maDiags = np.array([mainDiag, ludiag, ludiag])
#ma = spdiags(maDiags, diags, nPoint, nPoint).toarray()
#print("\ntridiagonal matrix ma:\n", ma)

# eigenvalues in ascending (aufsteigend) order, the normalized eigenvector
# corresponding to the eigenvalue eiva[i] is the column eive[:,i]
eiva, eive = eigh_tridiagonal(mainDiag, ludiag, select='a', select_range=None)
print("\neigenvalues\n", eiva, "\n\neigenvectors\n", eive)

# calculate squared probability density:
probDensity = np.square(np.absolute(eive))

# calculate normalized wavefunctions:
normFactors = np.zeros(nPoint)
eiveNorm = np.zeros((nPoint, nPoint))
probDensityNorm = np.zeros((nPoint, nPoint))
for ii in range(0, nPoint):
    normFactors[ii] = np.trapz(probDensity[:,ii], xinterp)
print("\nnormFactors\n", normFactors)
for jj in range(0, nPoint):
    eiveNorm[:,jj] = eive[:,jj] / np.sqrt(normFactors[jj])
    probDensityNorm[:,jj] = probDensity[:,jj] / normFactors[jj]
print("\nnormalized density should be 1\n", np.trapz(probDensityNorm[:,1], xinterp))
print("\nnormalized integrals should be 1\n", np.trapz(eiveNorm[:,2], xinterp))

# create wavefunction-points:
plt.plot(xinterp, yinterp, "b-", label="interp. Potential")
plt.plot(xPot, yPot, "ro", label="given Points")
plt.plot(xinterp, eiveNorm[:,2], "k--", label="normalized wave1")
#plt.plot(xinterp, eiveNorm[:,1], "k--", label="normalized wave2")
plt.plot(xinterp, probDensityNorm[:,2], "g--", label="normalized probDensity")
plt.xlim([xMin, xMax])
plt.ylim([-2,2])
plt.xlabel("$x$ [Bohr]")
plt.ylabel("Energy $E$ [Hartree]")
plt.title("Comparisation")
plt.legend()
plt.show()