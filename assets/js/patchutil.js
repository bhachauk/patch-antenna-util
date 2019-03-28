// All values are in mm.
// init area
var axis_i = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2]
var axis_j = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3]
var axis_k = [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6]

var ct= 0.03556;
var options = { "Copper" : "rgb(255, 140, 0)",
                "Silver" : "rgb(192, 192, 192)",
                "Gold" : "rgb(255, 215, 0)" }
var current_conductor = Object.keys(options)[0]

cavitycolor = getArray ('rgb(0, 100, 0)', 6)
coppercolor = getArray (options.Copper, 6)

cavityface = getFacecolor(cavitycolor);
copperface = getFacecolor(coppercolor);


var pw, pl, fl, fw, ch, tl, tfp, bfp, gl, gw, ele_l, effl, effd, dispMap;


function getConductorElement(){
    return document.getElementById("selectConductor");
}

function beReady(){
    var select = getConductorElement()
    if(select.value)
        return
    while(select.firstChild){
        select.removeChild(select.firstChild);
    }
    Object.keys(options).forEach(function(key){
        var el = document.createElement("option");
        el.textContent = key;
        el.value = key;
        select.add(el);
    })
}

function save(){
    var val = getConductorElement().value;
    var new_ct = parseFloat (document.getElementById('conductor_thickness_in').value)

    if(new_ct)
        ct = new_ct

    current_conductor = val;
    coppercolor = getArray(options[val], 6)
    copperface = getFacecolor(coppercolor);
    document.getElementById("info_btn").click();
}

function showInfo()
{
    setValue("current_conductor", "Current Conductor     :    " + current_conductor);
    setValue("conductor_thickness", "Conductor Thickness   :    " + ct + "   mm");
}


function updateParams() {
    tl = pl + fl
    tfp = (pw / 2) + (fw / 2)
    bfp = (pw / 2) - (fw / 2)
    dispMap ={
        "pw"   : ["Patch Width           ", pw],
        "pl"   : ["Patch Length          ", pl],
        "fw"   : ["Feeder Width          ", fw],
        "fl"   : ["Feeder Length         ", fl],
        "effl" : ["Effective Length      ", effl],
        "ele_l": ["Electrical Length     ", ele_l],
        "effd" : ["Effective Dielectric  ", effd],
        "gl"   : ["Ground Length         ", gl],
        "gw"   : ["Ground Width          ", gw]
    }
    for (key in dispMap){
        var val = dispMap[key]
        setValue( key, val.join([separator = ':  ']))
    }
}


// Patch calculations

function setValue(id, val)
{
    document.getElementById(id).innerHTML = val;
    document.getElementById(id).style.display= 'block';
}

function calculate(){
    var freq =  parseFloat (document.getElementById('f_in').value);
    var diel =  parseFloat (document.getElementById('d_in').value);
    var cavity =  parseFloat (document.getElementById('c_in').value);

    console.log("Got the frequency           : ",  freq)
    console.log("Got the dielectric constant : ",  diel)
    console.log("Got the frequency            : ",  cavity)

     if (! (freq && diel && cavity)){
        //simFailure();
        return false;
    }

    //all value in meter
    var ls = 299792458
    var f = freq * Math.pow(10, 9);
    var h = cavity * Math.pow(10, -3);

    // patch width
    var w = ( ls/(2*f)) * ( Math.sqrt (2/(diel + 1)) )
    console.log('width: ', w)

    // effective d
    effd = ((diel + 1)/2)*((diel-1)/2)*Math.pow((1+((12*h)/w)), -0.5)
    effl = (ls/(2*f))*(Math.sqrt(1/effd))

    var r = w/ h

    var dell = (0.412*h)*(( effd + 0.3)/(effd-0.264))*((r+0.258)/(r+0.8))

    var feed_l = (ls/(4*f))*(Math.sqrt(1/(effd)))

    var feed_w = (2*w)/5

    var patch_l = effl - (2 * dell)

    console.log('Ans: ', w, patch_l, feed_w, feed_l)

    ele_l = (360*feed_l)/ (ls/f)

    // Updating
    gl = (patch_l + feed_l + (6 * h)) * Math.pow(10, 3)
    gw = (w + feed_w + (6 * h))  * Math.pow(10, 3)
    ch = cavity
    pl = patch_l * Math.pow(10, 3)
    pw = w * Math.pow(10, 3)
    fl = feed_l * Math.pow(10, 3)
    fw = feed_w * Math.pow(10, 3)
    updateParams()
    //simSuccess()
    return true;
}

function simFailure()
{
    document.getElementById('simres').innerHTML = 'Some values Not entered ...';
    document.getElementById('simres').style.color = 'Red';
}

function simSuccess()
{
     document.getElementById('simres').innerHTML = 'Success';
     document.getElementById('simres').style.color = 'Green';
}

function getArray (v, n) {
    var arr = [];
    for (var i = 0; i < n; ++i) {
        arr[i] = v;
    }
    return arr;
}


function getFacecolor(facecolor){

    facecolor2 = new Array(facecolor.length * 2);
    facecolor.forEach(function(x, i) {
	    facecolor2[i * 2 + 1] = facecolor2[i * 2] = x;
    });
    return facecolor2;
}

function putChart() {

    var isvalid = calculate()
    if (!isvalid){
        return
    }

    var type = "mesh3d"
    var init =  { i: axis_i, j: axis_j, k: axis_k, type: "mesh3d"}
    var copperInit = Object.assign({}, init, {facecolor: copperface})

    var groundPts = getGroundPoints(tl, pw, ct)
    var ground = Object.assign({}, { x: groundPts[0], y: groundPts[1], z: groundPts[2], name: 'ground' }, copperInit)

    var cavityPts = getCavityPoints (tl, pw, ch, ct)
    var cavity = Object.assign({ x: cavityPts[0], y: cavityPts[1], z: cavityPts[2], name: 'cavity' }, init, {facecolor: cavityface})

    var patchPts = getPatchPoints (pl, pw, ch, ct)
    var patch = Object.assign({ x: patchPts[0], y: patchPts[1], z: patchPts[2], name: 'patch' }, copperInit)

    var feederPts = getFeederPoints (tl, pl, tfp, bfp, ch, ct)
    var feeder = Object.assign({ x: feederPts[0], y: feederPts[1], z: feederPts[2], name: 'feeder' }, copperInit)

    var data = [ ground, cavity, patch, feeder];

    var plt_r = tl * 1.2
    var layout = {

    plot_bgcolor:"#444",
    paper_bgcolor:"#444",
	scene: {
            xaxis:{title: 'Length'},
            yaxis:{title: 'Width'},
            zaxis:{title: 'Height'},

    xaxis: {
            nticks: 10,
            range: [-plt_r, 2 * plt_r]
            },
    yaxis: {
            nticks: 10,
            range: [-plt_r, 2 * plt_r]
            },
    zaxis: {
            nticks: 10,
            range: [-plt_r, plt_r]
            }
        }
    }
    Plotly.newPlot('plotd', data, layout, {});
    document.getElementById('param_btn').click();
}


function getGroundPoints (tl, pw, ct)
{
    var x = [0, 0, tl, tl, 0, 0, tl, tl]
    var y = [0, pw, pw, 0, 0, pw, pw, 0]
    var z = [0, 0, 0, 0, ct, ct, ct, ct]
    return [x,y,z]
}

function getCavityPoints (tl, pw, ch, ct)
{
    var x = [0, 0, tl, tl, 0, 0, tl, tl]
    var y = [0, pw, pw, 0, 0, pw, pw, 0]
    var z = [ct, ct, ct, ct, ch, ch, ch, ch]
    return [x,y,z]
}

function getPatchPoints (pl, pw, ch, ct)
{
    var x = [0, 0, pl, pl, 0, 0, pl, pl]
    var y = [0, pw, pw, 0, 0, pw, pw, 0]
    var z = [ch, ch, ch, ch, (ch + ct), (ch + ct), (ch + ct), (ch + ct)]
    return [x,y,z]
}

function getFeederPoints (tl, pl, tfp, bfp,ch, ct)
{
    var x = [ pl, pl, tl, tl, pl, pl, tl, tl]
    var y = [ bfp, tfp, tfp, bfp, bfp, tfp, tfp, bfp]
    var z = [ ch, ch, ch, ch, (ch + ct), (ch + ct), (ch + ct), (ch + ct)]
    return [x,y,z]
}

function toggle(evt, id) {
    makeDisplay ('show', 'none')
    x = document.getElementsByClassName('currentBtn');
    for (i = 0; i < x.length; i++) {
        x[i].className = x[i].className.replace("currentBtn", "cmnBtn");
    }
    var t = evt.currentTarget
    t.className = t.className.replace("cmnBtn", "currentBtn");
    makeDisplayId(id, 'block');
}

function makeDisplay (className, disType)
{
    var list = document.getElementsByClassName(className);
     for (i = 0; i < list.length; i++) {
       list[i].style.display = disType;
     }
}

function makeDisplayId (id, disType)
{
    document.getElementById (id).style.display = disType;
}
