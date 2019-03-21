import math
import matplotlib.pyplot as plt
import numpy as np
from math import cos, sin, sqrt, atan2, acos, pi, log10
import plotly
from plotly.offline import init_notebook_mode
import plotly.graph_objs as go
plotly.offline.init_notebook_mode(connected=True)
import scipy.integrate


def sph2cart1(r, th, phi):
  x = r * cos(phi) * sin(th)
  y = r * sin(phi) * sin(th)
  z = r * cos(th)
 
  return x, y, z
   
def cart2sph1(x, y, z):
  r = sqrt(x**2 + y**2 + z**2) + 1e-15
  th = acos(z / r)
  phi = atan2(y, x)
 
  return r, th, phi


def PatchFunction(thetaInDeg, phiInDeg, Freq, W, L, h, Er):
    """
    Taken from Design_patchr
    Calculates total E-field pattern for patch as a function of theta and phi
    Patch is assumed to be resonating in the (TMx 010) mode.
    E-field is parallel to x-axis
    W......Width of patch (m)
    L......Length of patch (m)
    h......Substrate thickness (m)
    Er.....Dielectric constant of substrate
    Refrence C.A. Balanis 2nd Edition Page 745
    """
    lamba = 3e8 / Freq

    theta_in = math.radians(thetaInDeg)
    phi_in = math.radians(phiInDeg)

    ko = 2 * math.pi / lamba

    xff, yff, zff = sph2cart1(999, theta_in, phi_in)                            # Rotate coords 90 deg about x-axis to match array_utils coord system with coord system used in the model.
    xffd = zff
    yffd = xff
    zffd = yff
    r, thp, php = cart2sph1(xffd, yffd, zffd)
    phi = php
    theta = thp

    if theta == 0:
        theta = 1e-9                                                              # Trap potential division by zero warning

    if phi == 0:
        phi = 1e-9

    Ereff = ((Er + 1) / 2) + ((Er - 1) / 2) * (1 + 12 * (h / W)) ** -0.5        # Calculate effictive dielectric constant for microstrip line of width W on dielectric material of constant Er

    F1 = (Ereff + 0.3) * (W / h + 0.264)                                        # Calculate increase length dL of patch length L due to fringing fields at each end, giving total effective length Leff = L + 2*dL
    F2 = (Ereff - 0.258) * (W / h + 0.8)
    dL = h * 0.412 * (F1 / F2)

    Leff = L + 2 * dL

    Weff = W                                                                    # Calculate effective width Weff for patch, uses standard Er value.
    heff = h * sqrt(Er)

    # Patch pattern function of theta and phi, note the theta and phi for the function are defined differently to theta_in and phi_in
    Numtr2 = sin(ko * heff * cos(phi) / 2)
    Demtr2 = (ko * heff * cos(phi)) / 2
    Fphi = (Numtr2 / Demtr2) * cos((ko * Leff / 2) * sin(phi))

    Numtr1 = sin((ko * heff / 2) * sin(theta))
    Demtr1 = ((ko * heff / 2) * sin(theta))
    Numtr1a = sin((ko * Weff / 2) * cos(theta))
    Demtr1a = ((ko * Weff / 2) * cos(theta))
    Ftheta = ((Numtr1 * Numtr1a) / (Demtr1 * Demtr1a)) * sin(theta)

    # Due to groundplane, function is only valid for theta values :   0 < theta < 90   for all phi
    # Modify pattern for theta values close to 90 to give smooth roll-off, standard model truncates H-plane at theta=90.
    # PatEdgeSF has value=1 except at theta close to 90 where it drops (proportional to 1/x^2) to 0

    rolloff_factor = 0.5                                                       # 1=sharp, 0=softer
    theta_in_deg = theta_in * 180 / math.pi                                          # theta_in in Deg
    F1 = 1 / (((rolloff_factor * (abs(theta_in_deg) - 90)) ** 2) + 0.001)       # intermediate calc
    PatEdgeSF = 1 / (F1 + 1)                                                    # Pattern scaling factor

    UNF = 1.0006                                                                # Unity normalisation factor for element pattern

    if theta_in <= math.pi / 2:
        Etot = Ftheta * Fphi * PatEdgeSF * UNF                                   # Total pattern by pattern multiplication
    else:
        Etot = 0

    return Etot

def GetPatchFields(PhiStart, PhiStop, ThetaStart, ThetaStop, Freq, W, L, h, Er):
    """"
    Calculates the E-field for range of thetaStart-thetaStop and phiStart-phiStop
    Returning a numpy array of form - fields[phiDeg][thetaDeg] = eField
    W......Width of patch (m)
    L......Length of patch (m)
    h......Substrate thickness (m)
    Er.....Dielectric constant of substrate
    """
    fields = np.ones((PhiStop, ThetaStop))                                      # Create initial array to hold e-fields for each position

    for phiDeg in range(PhiStart, PhiStop):
            for thetaDeg in range(ThetaStart, ThetaStop):                       # Iterate over all Phi/Theta combinations
                eField = PatchFunction(thetaDeg, phiDeg, Freq, W, L, h, Er)     # Calculate the field for current Phi, Theta
                fields[phiDeg][thetaDeg] = eField                               # Update array with e-field

    return fields


def PatchEHPlanePlot(Freq, W, L, h, Er, isLog=True):
    """
    Plot 2D plots showing E-field for E-plane (phi = 0) and the H-plane (phi = 90).
    """

    fields = GetPatchFields(0, 360, 0, 90, Freq, W, L, h, Er)                                                   # Calculate the field at each phi, theta

    Xtheta = np.linspace(0, 90, 90)                                                                             # Theta range array used for plotting

    if isLog:                                                                                                   # Can plot the log scale or normal
        plt.plot(Xtheta, 20 * np.log10(abs(fields[90, :])), label="H-plane (Phi=90)")                          # Log = 20 * log10(E-field)
        plt.plot(Xtheta, 20 * np.log10(abs(fields[0, :])), label="E-plane (Phi=0)")
        plt.ylabel('E-Field (dB)')
    else:
        plt.plot(Xtheta, fields[90, :], label="H-plane (Phi=90)")
        plt.plot(Xtheta, fields[0, :], label="E-plane (Phi=0)")
        plt.ylabel('E-Field')

    plt.xlabel('Theta (degs)')                                                                                 # Plot formatting
    plt.title("EH Plane - Theta ")
    plt.ylim(-40)
    plt.xlim((0, 90))

    start, end = plt.xlim()
    plt.xticks(np.arange(start, end, 5))
    plt.grid(b=True, which='major')
    plt.legend()
    plt.show()                                                                                                  # Show plot

    return fields                                                                                               # Return the calculated fields


def SurfacePlot(Fields, Freq, W, L, h, Er):
    """Plots 3D surface plot over given theta/phi range in Fields by calculating cartesian coordinate equivalent of spherical form."""

    print("Processing SurfacePlot...")

    phiSize = Fields.shape[0]                                                                                   # Finds the phi & theta range
    thetaSize = Fields.shape[1]

    X = np.ones((phiSize, thetaSize))                                                                           # Prepare arrays to hold the cartesian coordinate data.
    Y = np.ones((phiSize, thetaSize))
    Z = np.ones((phiSize, thetaSize))

    for phi in range(phiSize):                                                                                  # Iterate over all phi/theta range
        for theta in range(thetaSize):
            e = Fields[phi][theta]

            xe, ye, ze = sph2cart1(e, math.radians(theta), math.radians(phi))                                   # Calculate cartesian coordinates

            X[phi, theta] = xe                                                                                  # Store cartesian coordinates
            Y[phi, theta] = ye
            Z[phi, theta] = ze
    print X                                                                        # Plot surface
    # surface = go.Surface(x=X, y=Y, z=Z)
    # data = [surface]

    # layout = go.Layout(
    #     title='Surface Plot of EH Plane',
    #     scene=dict(
    #         xaxis=dict(
    #             gridcolor='rgb(255, 255, 255)',
    #             zerolinecolor='rgb(255, 255, 255)',
    #             showbackground=True,
    #             backgroundcolor='rgb(230, 230,230)'
    #         ),
    #         yaxis=dict(
    #             gridcolor='rgb(255, 255, 255)',
    #             zerolinecolor='rgb(255, 255, 255)',
    #             showbackground=True,
    #             backgroundcolor='rgb(230, 230,230)'
    #         ),
    #         zaxis=dict(
    #             gridcolor='rgb(255, 255, 255)',
    #             zerolinecolor='rgb(255, 255, 255)',
    #             showbackground=True,
    #             backgroundcolor='rgb(230, 230,230)'
    #         )
    #     )
    # )

    # fig = go.Figure(data=data, layout=layout)
    # plotly.offline.plot(fig)


def DesignPatch(Er, h, Freq):
    
    Eo = 8.854185e-12

    lambd = 3e8 / Freq
    lambdag = lambd / sqrt(Er)

    W = (3e8 / (2 * Freq)) * sqrt(2 / (Er + 1))

    temp = 1 + 12*(h/W)

    Ereff = ((Er + 1) / 2) + ((Er - 1) / 2) * temp ** -0.5                              # Calculate effictive dielectric constant for microstrip line of width W on dielectric material of constant Er

    F1 = (Ereff + 0.3) * (W / h + 0.264)                                                 # Calculate increase length dL of patch length L due to fringing fields at each end, giving actual length L = Lambda/2 - 2*dL
    F2 = (Ereff - 0.258) * (W / h + 0.8)
    dL = h * 0.412 * (F1 / F2)

    lambdag = lambd / sqrt(Ereff)
    L = (lambdag / 2) - 2 * dL

    print('Rectangular Microstrip Patch Design')
    print("Frequency: " + str(Freq))
    print("Dielec Const, Er : " + str(Er))
    print("Patch Width,  W: " + str(W) + "m")
    print("Patch Length,  L: " + str(L) + "m")
    print("Patch Height,  h: " + str(h) + "m")

    return W, L


def S_i(a):

    temp=scipy.integrate.quad(lambda x:sin(x)/x,0,a)

    return temp[0]

def J0(s):

    temp=scipy.integrate.quad(lambda x:cos(s*sin(x)),0,pi)

    temp=(1/pi)*temp[0]

    return temp

## Getting Conductance params

def getK0 (f):

    lamda_0 = v/f
    k0 = (2*pi)/lamda_0

    return k0


def getG1 (W, f):

    k0 = getK0 (f)
    X = k0 * W
    I1 = -2 + cos(X) + X*S_i(X) + sin(X)/X
    G1 = I1 / ( 120 * pi**2 )

    return G1

def getG12 (W, k0, L):

    temp=scipy.integrate.quad(lambda x:(((sin(k0*W*cos(x)/2)/cos(x))**2)*J0(k0*L*sin(x))*sin(x)**3),0,pi)
    G12=(1/(120*pi**2))*temp[0]
    return G12


def getGs(f, W, L):
    G1 = getG1(W, f)
    k0 = getK0(f)
    G12 = getG12(W, k0, L)
    return G1, G12

#-----------------------------------------#

def inputImpedance (f, W, L, h, epsilon):

    global v

    k0 = getK0 (f)
    G1, G12 = getGs(f, W, L)

    Rin=1/(2*(G1+G12))
    print "Input Impedance:",Rin,"ohms"
    return Rin


def insetFeedPosition(Rin, L):
    
    R=50.0
    
    y0=(L/pi)*(math.acos(sqrt(R/Rin)))
    
    return y0

def getDirectivity(G1, G12, W, f, I1, I2):

    global v

    lamda_0 = v/f

    g_12=G12/G1
    D_AF=2/(1+g_12)
    D0=((2*pi*W)/lamda_0)**2*(1/I1)
    D2=D0*D_AF

    DIR_1 = 10*log10(D2)

    I2=3.59801
    D_2=((2*pi*W)/lamda_0)**2*(pi/I2)
    DIR_2 = 10*log10(D_2)

    return DIR_1, DIR_2


v = 3 * 10 ** 8