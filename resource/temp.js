var f = 2.4 * Math.pow(10, 9)




function GetPatchFields(PhiStart, PhiStop, ThetaStart, ThetaStop, Freq, W, L, h, Er)

{
    fields = np.ones((PhiStop, ThetaStop))

    for phiDeg in range(PhiStart, PhiStop):
            for thetaDeg in range(ThetaStart, ThetaStop):
                eField = PatchFunction(thetaDeg, phiDeg, Freq, W, L, h, Er)
                fields[phiDeg][thetaDeg] = eField
    return fields
}

function PatchEHPlanePlot(Freq, W, L, h, Er){

    console.log(Freq, W, L, h, Er)
    fields = GetPatchFields(0, 360, 0, 90, Freq, W, L, h, Er)
}