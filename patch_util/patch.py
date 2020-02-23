import math
import matplotlib.pyplot as plt
import numpy as np
from math import cos, sin, sqrt, atan2, acos, pi, log10
import plotly
from plotly.offline import iplot
import scipy.integrate
import plotly.graph_objs as go
plotly.offline.init_notebook_mode(connected=True)

v = 3 * 10 ** 8


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


def patch_function(theta_in_deg, phi_in_deg, freq, w, l, h, er):
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
    lambda_ = 3e8 / freq
    theta_in = math.radians(theta_in_deg)
    phi_in = math.radians(phi_in_deg)

    ko = 2 * math.pi / lambda_

    xff, yff, zff = sph2cart1(999, theta_in, phi_in)                            # Rotate coords 90 deg about x-axis to match array_utils coord system with coord system used in the model.
    xffd = zff
    yffd = xff
    zffd = yff
    r, thp, php = cart2sph1(xffd, yffd, zffd)
    phi = php
    theta = thp

    if theta == 0:
        # Trap potential division by zero warning
        theta = 1e-9

    if phi == 0:
        phi = 1e-9

    # Calculate effective dielectric constant for micro_strip line of width W on dielectric material of constant Er
    e_ref = ((er + 1) / 2) + ((er - 1) / 2) * (1 + 12 * (h / w)) ** -0.5

    # Calculate increase length dL of patch length L due to fringing fields at each end,
    # giving total effective length Leff = L + 2*dL

    f1 = (e_ref + 0.3) * (w / h + 0.264)
    f2 = (e_ref - 0.258) * (w / h + 0.8)
    d_l = h * 0.412 * (f1 / f2)

    l_eff = l + 2 * d_l

    # Calculate effective width Weff for patch, uses standard Er value.
    w_eff = w
    h_eff = h * sqrt(er)

    # Patch pattern function of theta and phi,
    # Note the theta and phi for the function are defined differently to theta_in and phi_in
    num_tr_2 = sin(ko * h_eff * cos(phi) / 2)
    dem_tr_2 = (ko * h_eff * cos(phi)) / 2
    f_phi = (num_tr_2 / dem_tr_2) * cos((ko * l_eff / 2) * sin(phi))

    num_tr_1 = sin((ko * h_eff / 2) * sin(theta))
    dem_tr_1 = ((ko * h_eff / 2) * sin(theta))
    num_tr_1a = sin((ko * w_eff / 2) * cos(theta))
    dem_tr_1a = ((ko * w_eff / 2) * cos(theta))
    f_theta = ((num_tr_1 * num_tr_1a) / (dem_tr_1 * dem_tr_1a)) * sin(theta)

    # Due to groundplane, function is only valid for theta values :   0 < theta < 90   for all phi
    # Modify pattern for theta values close to 90 to give smooth roll-off, standard model truncates H-plane at theta=90.
    # PatEdgeSF has value=1 except at theta close to 90 where it drops (proportional to 1/x^2) to 0

    # 1=sharp, 0=softer
    roll_off_factor = 0.5
    # theta_in in Deg
    theta_in_deg = theta_in * 180 / math.pi
    # intermediate calc
    f1 = 1 / (((roll_off_factor * (abs(theta_in_deg) - 90)) ** 2) + 0.001)
    # Pattern scaling factor
    pat_edge_sf = 1 / (f1 + 1)
    # Unity normalisation factor for element pattern
    UNF = 1.0006

    # Total pattern by pattern multiplication
    if theta_in <= math.pi / 2:
        e_tot = f_theta * f_phi * pat_edge_sf * UNF
    else:
        e_tot = 0

    return e_tot


def get_patch_fields(phi_start, phi_stop, theta_start, theta_stop, freq, w, l, h, er):
    """"
    Calculates the E-field for range of thetaStart-thetaStop and phiStart-phiStop
    Returning a numpy array of form - fields[phiDeg][thetaDeg] = eField
    W......Width of patch (m)
    L......Length of patch (m)
    h......Substrate thickness (m)
    Er.....Dielectric constant of substrate
    """
    # Create initial array to hold e-fields for each position
    fields = np.ones((phi_stop, theta_stop))
    # Iterate over all Phi/Theta combinations
    for phiDeg in range(phi_start, phi_stop):
        for thetaDeg in range(theta_start, theta_stop):
            # Calculate the field for current Phi, Theta
            eField = patch_function(thetaDeg, phiDeg, freq, w, l, h, er)
            # Update array with e-field
            fields[phiDeg][thetaDeg] = eField

    return fields


def patch_eh_plane_plot(freq, w, l, h, er, is_log=True):
    """
    Plot 2D plots showing E-field for E-plane (phi = 0) and the H-plane (phi = 90).
    """

    fields = get_patch_fields(0, 360, 0, 90, freq, w, l, h, er)

    Xtheta = np.linspace(0, 90, 90)

    if is_log:
        # Log = 20 * log10(E-field)# Can plot the log scale or normal
        plt.plot(Xtheta, 20 * np.log10(abs(fields[90, :])), label="H-plane (Phi=90)")
        plt.plot(Xtheta, 20 * np.log10(abs(fields[0, :])), label="E-plane (Phi=0)")
        plt.ylabel('E-Field (dB)')
    else:
        plt.plot(Xtheta, fields[90, :], label="H-plane (Phi=90)")
        plt.plot(Xtheta, fields[0, :], label="E-plane (Phi=0)")
        plt.ylabel('E-Field')

    plt.xlabel('Theta (degs)')
    plt.title("EH Plane - Theta ")
    plt.ylim(-40)
    plt.xlim((0, 90))

    start, end = plt.xlim()
    plt.xticks(np.arange(start, end, 5))
    plt.grid(b=True, which='major')
    plt.legend()
    plt.show()
    return fields


def surface_plot(fields, is_note_book=False):
    """Plots 3D surface plot over given theta/phi range in Fields by calculating cartesian
    coordinate equivalent of spherical form."""

    print("Processing SurfacePlot...")
    # Finds the phi & theta range
    phiSize = fields.shape[0]
    thetaSize = fields.shape[1]
    # Prepare arrays to hold the cartesian coordinate data.
    X = np.ones((phiSize, thetaSize))
    Y = np.ones((phiSize, thetaSize))
    Z = np.ones((phiSize, thetaSize))

    for phi in range(phiSize):
        for theta in range(thetaSize):
            e = fields[phi][theta]

            xe, ye, ze = sph2cart1(e, math.radians(theta), math.radians(phi))

            X[phi, theta] = xe
            Y[phi, theta] = ye
            Z[phi, theta] = ze
    surface = go.Surface(x=X, y=Y, z=Z)
    data = [surface]

    layout = go.Layout(
        title='Surface Plot of EH Plane',
        scene=dict(
            xaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            yaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            zaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            )
        )
    )

    fig = go.Figure(data=data, layout=layout)
    if is_note_book:
        iplot(fig)
    else:
        plotly.offline.plot(fig)


def design_patch(er, h, freq):

    lambda_ = 3e8 / freq
    w = (3e8 / (2 * freq)) * sqrt(2 / (er + 1))
    temp = 1 + 12*(h/w)

    e_ref = ((er + 1) / 2) + ((er - 1) / 2) * temp ** -0.5

    f1 = (e_ref + 0.3) * (w / h + 0.264)
    f2 = (e_ref - 0.258) * (w / h + 0.8)
    d_l = h * 0.412 * (f1 / f2)

    lambda_g = lambda_ / sqrt(e_ref)
    L = (lambda_g / 2) - 2 * d_l

    print('Rectangular Microstrip Patch Design')
    print("Frequency: " + str(freq))
    print("Dielec Const, Er : " + str(er))
    print("Patch Width,  W: " + str(w) + "m")
    print("Patch Length,  L: " + str(L) + "m")
    print("Patch Height,  h: " + str(h) + "m")
    return w, L


def S_i(a):
    temp = scipy.integrate.quad(lambda x:sin(x)/x,0,a)
    return temp[0]


def J0(s):
    temp = scipy.integrate.quad(lambda x:cos(s*sin(x)),0,pi)
    temp = (1/pi)*temp[0]
    return temp


def get_k(f):
    lamda_0 = v/f
    k0 = (2*pi)/lamda_0
    return k0


def getG1 (W, f):
    k0 = get_k (f)
    X = k0 * W
    I1 = -2 + cos(X) + X*S_i(X) + sin(X)/X
    G1 = I1 / ( 120 * pi**2 )
    return G1


def getG12 (W, k0, L):
    temp = scipy.integrate.quad(lambda x: (((sin(k0*W*cos(x)/2)/cos(x))**2)*J0(k0*L*sin(x))*sin(x)**3), 0, pi)
    G12 = (1/(120*pi**2))*temp[0]
    return G12


def getGs(f, W, L):
    G1 = getG1(W, f)
    k0 = get_k(f)
    G12 = getG12(W, k0, L)
    return G1, G12


def input_impedance (f, W, L):

    global v
    k0 = get_k (f)
    G1, G12 = getGs(f, W, L)
    Rin = 1/(2*(G1+G12))
    print("Input Impedance:", Rin, "ohms")
    return Rin


def inset_feed_position(Rin, L):
    
    R = 50.0
    y0 = (L/pi)*(math.acos(sqrt(R/Rin)))
    return y0


def get_directivity(G1, G12, W, f, I1, I2):
    global v
    lamda_0 = v/f
    g_12 = G12/G1
    D_AF = 2/(1+g_12)
    D0 = ((2*pi*W)/lamda_0)**2*(1/I1)
    D2 = D0 * D_AF
    DIR_1 = 10*log10(D2)
    D_2 = ((2*pi*W)/lamda_0) ** 2 * (pi/I2)
    DIR_2 = 10 * log10(D_2)
    return DIR_1, DIR_2
