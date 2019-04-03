var f = 2.4 * Math.pow(10, 9)

var toMM = Math.pow(10, -3)

var pw = 23.700674815509178 * toMM
var pl = 21.587764623402762 * toMM
var fw = 9.480269926203672 * toMM
var fl = 12.760735864844596 * toMM
var effl = 0.02552147172968919 * toMM
var ele_l = 61.29393570460074
var ed = 2.156006183008272 * toMM
var gl = 58.34850048824735 * toMM
var gw = 57.18094474171285 * toMM
var h = 1.6 * toMM
var er = 4.4

var Fields = GetPatchFields(0, 360, 0, 90)

//console.log('fields', Fields)

phiSize = 360
thetaSize = 90

X = Array(phiSize).fill(Array(thetaSize).fill(1));
Y = Array(phiSize).fill(Array(thetaSize).fill(1));
Z = Array(phiSize).fill(Array(thetaSize).fill(1));

for (var phi_temp =0; phi_temp < phiSize; phi_temp++)
{
    for (var theta_temp =0; theta_temp < thetaSize; theta_temp++)
    {
        var e_temp = Fields[phi_temp][theta_temp]
        var t_temp = toRadians(theta_temp)
        var p_temp = toRadians(phi_temp)
        var ax = sph2cart1 (e_temp, t_temp, p_temp)
        X[phi_temp][theta_temp] = ax[0]
        Y[phi_temp][theta_temp] = ax[1]
        Z[phi_temp][theta_temp] = ax[2]
    }
}

var layout = {
    title: 'Surface Plot of EH Plane',
    scene: {
        xaxis: {
            gridcolor: 'rgb(255, 255, 255)',
            zerolinecolor: 'rgb(255, 255, 255)',
            showbackground: true,
            backgroundcolor: 'rgb(230, 230,230)'
        },
        yaxis: {
            gridcolor: 'rgb(255, 255, 255)',
            zerolinecolor: 'rgb(255, 255, 255)',
            showbackground: true,
            backgroundcolor: 'rgb(230, 230,230)'
        },
        zaxis: {
            gridcolor: 'rgb(255, 255, 255)',
            zerolinecolor: 'rgb(255, 255, 255)',
            showbackground: true,
            backgroundcolor: 'rgb(230, 230,230)'
        }
    }
}

var data = [{
        x: X,
        y: Y,
        z: Z,
        type: 'surface'
}]

function showPlot(){
    document.getElementById('plotd').display='block';
    Plotly.newPlot('plotd', data=data, {}, {});
}

function GetPatchFields(PhiStart, PhiStop, ThetaStart, ThetaStop)

{
    fields = Array(PhiStop).fill(Array(ThetaStop).fill(1));

    for (var phiDeg = PhiStart; phiDeg < PhiStop; phiDeg++) {
        for (var thetaDeg = ThetaStart; thetaDeg < ThetaStop; thetaDeg ++){
            var eField = PatchFunction (thetaDeg, phiDeg)
            fields[phiDeg][thetaDeg] = eField
        }
    }

    return fields
}

function toRadians (angle) {
  return angle * (Math.PI / 180);
}

function sph2cart1(r, th, phi){
  x = r * Math.cos(phi) * Math.sin(th)
  y = r * Math.sin(phi) * Math.sin(th)
  z = r * Math.cos(th)
  return [x, y, z]
}

function cart2sph1(x, y, z){
  r = Math.sqrt( Math.pow(x,2) + Math.pow(y,2) + Math.pow(z,2)) + Math.pow(10, -15)
  th = Math.acos(z / r)
  phi = Math.atan2(y, x)
  return [r, th, phi]
}

function PatchFunction (thetaInDeg, phiInDeg){

    var lambda = 3 * Math.pow(10, 8) / f;

    var theta_in = toRadians(thetaInDeg)
    var phi_in = toRadians(phiInDeg)

    var ko = 2 * Math.PI / lambda

    ff = sph2cart1(999, theta_in, phi_in)
    out = cart2sph1(ff[2], ff[0], ff[1])
    r = out[0]; theta = out[1]; phi = out[2];

    if (theta == 0)
       theta = Math.pow(10, -9)

    if (phi == 0)
        phi = Math.pow(10, -9)

    Ereff = ((er + 1) / 2) + ((er - 1) / 2) * (1 + 12 * (h / pw))
    Ereff = Math.pow(Ereff, -0.5)

    F1 = (Ereff + 0.3) * (pw / h + 0.264)
    F2 = (Ereff - 0.258) * (pw / h + 0.8)
    dL = h * 0.412 * (F1 / F2)

    Leff = pl + 2 * dL

    Weff = pw
    heff = h * Math.sqrt(er)

    Numtr2 = Math.sin(ko * heff * Math.cos(phi) / 2)
    Demtr2 = (ko * heff * Math.cos(phi)) / 2
    Fphi = (Numtr2 / Demtr2) * Math.cos((ko * Leff / 2) * Math.sin(phi))

    Numtr1 = Math.sin((ko * heff / 2) * Math.sin(theta))
    Demtr1 = ((ko * heff / 2) * Math.sin(theta))
    Numtr1a = Math.sin((ko * Weff / 2) * Math.cos(theta))
    Demtr1a = ((ko * Weff / 2) * Math.cos(theta))

    Ftheta = ((Numtr1 * Numtr1a) / (Demtr1 * Demtr1a)) * Math.sin(theta)

    rolloff_factor = 0.5
    theta_in_deg = theta_in * 180 / Math.PI
    F1 = 1 / (((rolloff_factor * (Math.abs(theta_in_deg) - 90)) ** 2) + 0.001)
    PatEdgeSF = 1 / (F1 + 1)

    UNF = 1.0006

    if (theta_in <= Math.PI / 2)
        Etot = Ftheta * Fphi * PatEdgeSF * UNF
    else
        Etot = 0

    return Etot
}