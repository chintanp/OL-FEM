// Variable declaration
// Some of these are matrices, so they will need to treated separately

var sylvester = require('sylvester'),
	Matrix = sylvester.Matrix,
	Vector = sylvester.Vector;

var mathjs = require('mathjs'),
    math = mathjs();


var nnd, nel, nne, nodof, eldof, n, ngpb, ngps;
var dim;
var W, MX, MY, MXY, QX, QY;

var geom = math.matrix();
var connec = math.matrix();
var deeb = math.matrix();
var dees = math.matrix();
var nf = math.matrix();
var load = math.matrix();
var sampb = math.matrix();
var samps = math.matrix();


// ******** These are the dimensions of the quarter plate, for total dimension of plate multiply each by two.******
var Lx = 6; 			// Test values, will update from UI later on
var Ly = 8;
var thick = 0.2;

var divx = 2;			// Divisions in X and Y
var divy = 2;

var etype_index = 2;    // Type of element, 1=4-noded or 2=8-noded
var end_index = 1;      // Type of end condition, 1=SS, 2=Fixed, 3=Beams
var loadtype_index = 1;     // Type of load - 1=Conc, 2=UDL

var X_origin = 0;
var Y_origin = 0;
var dim = 2;

// Functions - all right here, later need to be placed in separate files and referenced here

// This function generates the 4-noded bi-linear element
var Q4_mesh = function(len, wid, div_x, div_y, xo, yo) {

	var Length = len;
	var Width = wid;
	var NXE = div_x;
	var NYE = div_y;
	var X_origin = xo;
	var Y_origin = yo;
	var dhx = Length / NXE;
	var dhy = Width / NYE;

	this.geom = math.matrix();
    this.connec = math.matrix();
	this.nel = 0;
	this.nnd = 0;
	var k = 0;

	for( var i = 1; i <= NXE; i++ )
	{
		for( var j = 1; j <= NYE; j++ )
		{
			k = k + 1;

            n1 = j + (i-1) * (NYE+1);
            n2 = j + i * (NYE+1);
            n3 = n1 + 1;
            n4 = n2 + 1;

            this.geom.subset(math.index(n1-1, 0), (i-1)*dhx-X_origin);
            this.geom.subset(math.index(n1-1, 1), (j-1)*dhy-Y_origin);
            this.geom.subset(math.index(n2-1, 0), (i)*dhx-X_origin);
            this.geom.subset(math.index(n2-1, 1), (j-1)*dhy-Y_origin);
            this.geom.subset(math.index(n3-1, 0), (i-1)*dhx-X_origin);
            this.geom.subset(math.index(n3-1, 1), (j)*dhy-Y_origin);
            this.geom.subset(math.index(n4-1, 0), (i)*dhx-X_origin);
            this.geom.subset(math.index(n4-1, 1), (j)*dhy-Y_origin);

            this.nel = k;
            this.nnd = n4;

            this.connec.subset(math.index(this.nel-1, 0), n1);
            this.connec.subset(math.index(this.nel-1, 1), n2);
            this.connec.subset(math.index(this.nel-1, 2), n4);
            this.connec.subset(math.index(this.nel-1, 3), n3);

        }
	}

    return {nel: this.nel, nnd: this.nnd, geom: this.geom, connec: this.connec};
}

var Q8_mesh = function(len, wid, div_x, div_y, xo, yo) {

    var Length = len;
    var Width = wid;
    var NXE = div_x;
    var NYE = div_y;
    var X_origin = xo;
    var Y_origin = yo;
    var dhx = Length / NXE;
    var dhy = Width / NYE;

    this.geom = math.matrix();
    this.connec = math.matrix();
    this.nel = 0;
    this.nnd = 0;
    var k = 0;

    for( var i = 1; i <= NXE; i++ )
    {
        for( var j = 1; j <= NYE; j++ )
        {
            k = k + 1;

            n1 = (i-1) * (3*NYE + 2) + 2*j - 1;
            n2 = i * (3*NYE + 2) + j - NYE - 1;
            n3 = i * (3*NYE + 2) + 2*j - 1;
            n4 = n3 + 1;
            n5 = n3 + 2;
            n6 = n2 + 1;
            n7 = n1 + 2;
            n8 = n1 + 1;

            this.geom.subset(math.index(n1-1, 0), (i-1)*dhx-X_origin);
            this.geom.subset(math.index(n1-1, 1), (j-1)*dhy-Y_origin);
            this.geom.subset(math.index(n3-1, 0), (i)*dhx-X_origin);
            this.geom.subset(math.index(n3-1, 1), (j-1)*dhy-Y_origin);

            this.geom.subset(math.index(n2-1, 0), (this.geom.subset(math.index(n1-1, 0))+this.geom.subset(math.index(n3-1, 0)))/2);
            this.geom.subset(math.index(n2-1, 1), (this.geom.subset(math.index(n1-1, 1))+this.geom.subset(math.index(n3-1, 1)))/2);

            this.geom.subset(math.index(n5-1, 0), (i)*dhx-X_origin);
            this.geom.subset(math.index(n5-1, 1), (j)*dhy-Y_origin);

            this.geom.subset(math.index(n4-1, 0), (this.geom.subset(math.index(n3-1, 0))+this.geom.subset(math.index(n5-1, 0)))/2);
            this.geom.subset(math.index(n4-1, 1), (this.geom.subset(math.index(n3-1, 1))+this.geom.subset(math.index(n5-1, 1)))/2);

            this.geom.subset(math.index(n7-1, 0), (i-1)*dhx-X_origin);
            this.geom.subset(math.index(n7-1, 1), (j)*dhy-Y_origin);

            this.geom.subset(math.index(n6-1, 0), (this.geom.subset(math.index(n5-1, 0))+this.geom.subset(math.index(n7-1, 0)))/2);
            this.geom.subset(math.index(n6-1, 1), (this.geom.subset(math.index(n5-1, 1))+this.geom.subset(math.index(n7-1, 1)))/2);

            this.geom.subset(math.index(n8-1, 0), (this.geom.subset(math.index(n1-1, 0))+this.geom.subset(math.index(n7-1, 0)))/2);
            this.geom.subset(math.index(n8-1, 1), (this.geom.subset(math.index(n1-1, 1))+this.geom.subset(math.index(n7-1, 1)))/2);

            this.nel = k;
            this.nnd = n5;

            this.connec.subset(math.index(this.nel-1, 0), n1);
            this.connec.subset(math.index(this.nel-1, 1), n2);
            this.connec.subset(math.index(this.nel-1, 2), n3);
            this.connec.subset(math.index(this.nel-1, 3), n4);
            this.connec.subset(math.index(this.nel-1, 4), n5);
            this.connec.subset(math.index(this.nel-1, 5), n6);
            this.connec.subset(math.index(this.nel-1, 6), n7);
            this.connec.subset(math.index(this.nel-1, 7), n8);
        }
    }

    return {nel: this.nel, nnd: this.nnd, geom: this.geom, connec: this.connec};
}

var formdeeb = function(E, vu, thick) {

    var deeb = math.matrix();

    DR = E * Math.pow(thick, 3)/(12 * (1.0 - vu*vu));

    deeb.subset(math.index(0, 0), DR*1);
    deeb.subset(math.index(0, 1), DR*vu);
    deeb.subset(math.index(0, 2), DR*0);
    deeb.subset(math.index(1, 0), DR*vu);
    deeb.subset(math.index(1, 1), DR*1);
    deeb.subset(math.index(1, 2), DR*0);
    deeb.subset(math.index(2, 0), DR*0);
    deeb.subset(math.index(2, 1), DR*0);
    deeb.subset(math.index(2, 2), DR*(1.0 - vu)/2);

    return deeb;
}

var formdees = function(E, vu, thick) {

    var dees = math.matrix();

    G = E / (2 * (1.0+vu));

    dees.subset(math.index(0, 0), G*thick);
    dees.subset(math.index(0, 1), G*0);
    dees.subset(math.index(1, 0), G*0);
    dees.subset(math.index(1, 1), G*thick);

    return dees;
}

var gauss = function(ngp) {

    // This function returns the abscissa and weights of the Gauss points for ngp upto 4

    var samp = math.zeros(ngp,2);

    if(ngp == 1)
    {
        samp.subset(math.index(0,0), 0.0);
        samp.subset(math.index(0,1), 2.0);
    }
    else if (ngp == 2)
    {
        samp.subset(math.index(0,0), -1.0/Math.sqrt(3));
        samp.subset(math.index(0,1), 1.0);
        samp.subset(math.index(1,0), 1.0/Math.sqrt(3));
        samp.subset(math.index(1,1), 1.0);
    }
    else if (ngp == 3)
    {
        samp.subset(math.index(0,0), -0.2*Math.sqrt(15.0));
        samp.subset(math.index(0,1), 5.0/9.0);
        samp.subset(math.index(1,0), 0.0);
        samp.subset(math.index(1,1), 8.0/9.0);
        samp.subset(math.index(2,0), 0.2*Math.sqrt(15.0));
        samp.subset(math.index(2,1), 5.0/9.0);
    }
    else if (ngp == 4)
    {
        samp.subset(math.index(0,0), -0.861136311594053);
        samp.subset(math.index(0,1), 0.347854845137454);
        samp.subset(math.index(1,0), -0.339981043584856);
        samp.subset(math.index(1,1), 0.652145154862546);
        samp.subset(math.index(2,0), 0.339981043584856);
        samp.subset(math.index(2,1), 0.652145154862546);
        samp.subset(math.index(3,0), 0.861136311594053);
        samp.subset(math.index(3,1), 0.347854845137454);
    }
    return samp;

}

var platelem_q8 = function(input, i) {

    // This function returns the coordinates of the nodes of element i
    // and its steering vector
    this.coord = math.zeros(input.nne, input.dim);
    this.g = math.matrix();

    for(var k = 1; k <= input.nne; k++)
    {
        for(var j = 1; j <= input.dim; j++)
        {
            this.coord.subset(math.index(k-1,j-1), input.geom.subset(math.index(input.connec.subset(math.index(i-1,k-1))-1,j-1)));
        }
    }
    var l = 0;
    for(var k = 1; k <= input.nne; k++)
    {
        for(var j = 1; j <= input.nodof; j++)
        {
            l = l + 1;
            this.g.subset(math.index(l-1), input.nf.subset(math.index(input.connec.subset(math.index(i-1,k-1))-1, j-1)))
        }
    }
    return {coord: this.coord, g: this.g};
}

var fmlin = function(samp, ig, jg) {

    // This function returns the vector of shape function and its derivative
    // for 4-noded bi-linear quadrilateral element
    var xi = samp.subset(math.index(ig-1,0));
    var eta = samp.subset(math.index(jg-1,0));

    var fun = math.matrix();
    var der = math.matrix();

    fun.subset(math.index(0,0), 0.25*(1.0-xi-eta+xi*eta));
    fun.subset(math.index(1,0), 0.25*(1.0+xi-eta-xi*eta));
    fun.subset(math.index(2,0), 0.25*(1.0+xi+eta+xi*eta));
    fun.subset(math.index(3,0), 0.25*(1.0-xi+eta-xi*eta));

    der.subset(math.index(0,0), 0.25*(eta-1));
    der.subset(math.index(0,1), 0.25*(1-eta));
    der.subset(math.index(0,2), 0.25*(eta+1));
    der.subset(math.index(0,3), -0.25*(eta+1));
    der.subset(math.index(1,0), 0.25*(xi-1));
    der.subset(math.index(1,1), -0.25*(xi+1));
    der.subset(math.index(1,2), 0.25*(xi+1));
    der.subset(math.index(1,3), 0.25*(1-xi));

    return {der: this.der, fun: this.fun};
}

var fmquad = function(samp, ig, jg) {

    // This function returns the vector of the shape function and their derivatives
    // with respect to xi and eta at the Gauss points for an 8-noded quadrilateral
    var xi = samp.subset(math.index(ig-1,0));
    var eta = samp.subset(math.index(jg-1,0));
    var etam = (1.0 - eta);
    var etap = (1.0 + eta);
    var xim = (1.0 - xi);
    var xip = (1.0 + xi);

    this.der = math.matrix();
    this.fun = math.matrix();

    this.fun.subset(math.index(0), -0.25*xim*etam*(1.0 + xi + eta));
    this.fun.subset(math.index(1), 0.5*(1.0 - Math.pow(xi,2))*etam);
    this.fun.subset(math.index(2), -0.25*xip*etam*(1.0 - xi + eta));
    this.fun.subset(math.index(3), 0.5*xip*(1.0 - Math.pow(eta,2)));
    this.fun.subset(math.index(4), -0.25*xip*etap*(1.0 - xi - eta));
    this.fun.subset(math.index(5), 0.5*(1.0 - Math.pow(xi,2))*etap);
    this.fun.subset(math.index(6), -0.25*xim*etap*(1.0 + xi - eta));
    this.fun.subset(math.index(7), 0.5*xim*(1.0 - Math.pow(eta,2)));

    this.der.subset(math.index(0,0), 0.25*etam*(2.0*xi + eta));
    this.der.subset(math.index(0,1), -1.0*etam*xi);
    this.der.subset(math.index(0,2), 0.25*etam*(2.0*xi - eta));
    this.der.subset(math.index(0,3), 0.5*(1 - Math.pow(eta,2)));
    this.der.subset(math.index(0,4), 0.25*etap*(2.0*xi + eta));
    this.der.subset(math.index(0,5), -1.0*etap*xi);
    this.der.subset(math.index(0,6), 0.25*etap*(2.0*xi - eta));
    this.der.subset(math.index(0,7), -0.5*(1.0 - Math.pow(eta,2)));
    this.der.subset(math.index(1,0), 0.25*xim*(2.0*eta + xi));
    this.der.subset(math.index(1,1), -0.5*(1.0 - Math.pow(xi,2)));
    this.der.subset(math.index(1,2), -0.25*xip*(xi - 2.0*eta));
    this.der.subset(math.index(1,3), -1.0*xip*eta);
    this.der.subset(math.index(1,4), 0.25*xip*(xi + 2.0*eta));
    this.der.subset(math.index(1,5), 0.5*(1.0 - Math.pow(xi,2)));
    this.der.subset(math.index(1,6), -0.25*xim*(xi - 2.0*eta));
    this.der.subset(math.index(1,7), -1.0*xim*eta);

    return {der: this.der, fun: this.fun};

}

var formbeeb = function(deriv, nne, eldof) {

    // This function assembles the matrix [beeb] from the derivatives of the
    // shape functions in global coordinates for a thick plate element (bending action)

    var beeb = math.zeros(3, eldof);
    var k = 0;
    var j = 0;
    var x = 0;
    var y = 0;


    for(var m = 1; m <= nne; m++)
    {
        k = 3 * m;
        j = k - 1;
        x = deriv.subset(math.index(0,m-1));
        beeb.subset(math.index(0,j-1), x);
        beeb.subset(math.index(2, k-1), x);
        y = deriv.subset(math.index(1,m-1));
        beeb.subset(math.index(1,k-1), y);
        beeb.subset(math.index(2, j-1), y);
    }
    return beeb;

}

var formbees = function(deriv, fun, nne, eldof) {

    // This function assembles the matrix [bees] from the derivatives of the
    // shape functions in global coordinates for the shear action in plate element

    var bees = math.zeros(2, eldof);
    var k = 0;
    var j = 0;
    var x = 0;
    var y = 0;
    var i = 0;


    for(var m = 1; m <= nne; m++)
    {
        k = 3 * m;
        j = k - 1;
        i = k - 2;
        x = deriv.subset(math.index(0,m-1));
        y = deriv.subset(math.index(1,m-1));
        bees.subset(math.index(1,i-1), -1*x);
        bees.subset(math.index(0,i-1), -1*y);
        bees.subset(math.index(0,k-1), fun.subset(math.index(m-1)));
        bees.subset(math.index(1,j-1), fun.subset(math.index(m-1)));
    }
    return bees;

}
var form_KK = function(input, KK, kg, g) {

    // This function assembles the global stiffness matrix

    var KK = math.matrix();

    for(var i = 1; i <= input.eldof; i++)
    {
        var gi = g.subset(math.index(i-1));
        if(gi != 0)
        {
            for(var j = 1; j <= input.eldof; j++)
            {
                var gj = g.subset(math.index(j-1));
                if(gj != 0)
                {
                    KK.subset(math.index(gi-1, gj-1), KK.subset(math.index(gi-1, gj-1)) + kg.subset(math.index(i-1,j-1)));
                }
            }
        }
    }

    return KK;
}
// ******************************************************************************
// ******************************************************************************
if (etype_index == 1)           // 4-noded thick plate
{
    nne = 4;

    nel = Q4_mesh(Lx, Ly, divx, divy, X_origin, Y_origin).nel;
    nnd = Q4_mesh(Lx, Ly, divx, divy, X_origin, Y_origin).nnd;
    geom = Q4_mesh(Lx, Ly, divx, divy, X_origin, Y_origin).geom;
    connec = Q4_mesh(Lx, Ly, divx, divy, X_origin, Y_origin).connec;

    ngpb = 2;
    ngps = 1;
}
else if(etype_index == 2)           // 8-noded thick plate
{
    nne = 8;

    nel = Q8_mesh(Lx, Ly, divx, divy, X_origin, Y_origin).nel;
    nnd = Q8_mesh(Lx, Ly, divx, divy, X_origin, Y_origin).nnd;
    geom = Q8_mesh(Lx, Ly, divx, divy, X_origin, Y_origin).geom;
    connec = Q8_mesh(Lx, Ly, divx, divy, X_origin, Y_origin).connec;

    ngpb = 3;
    ngps = 2;
}

nodof = 3;              // Number of degrees of freedom per node
eldof = nne * nodof;    // number of degrees of freedom per element

density = 0.090318230000872 ;
E = 30.e+6;
vu = 0.3;

deeb = formdeeb(E, vu, thick);
dees = formdees(E, vu, thick);

// ********************************************
// Boundary Conditions

nf = math.ones(nnd, nodof);         // each row of nf corresponds to w_z, theta_x, and theta_y

if(end_index == 1)              // Simply supported Ends
{
    for (var i = 1; i <= nnd; i++)
    {
        if (geom.subset(math.index(i - 1, 0)) == 0)
        {
            nf.subset(math.index(i - 1, 0), 0);
            nf.subset(math.index(i - 1, 2), 0);
        }
        else if (geom.subset(math.index(i - 1, 1)) == 0)
        {
            nf.subset(math.index(i - 1, 0), 0);
            nf.subset(math.index(i - 1, 1), 0);
        }
        else if (geom.subset(math.index(i - 1, 0)) == Lx)
        {
            nf.subset(math.index(i - 1, 1), 0);
        }
        else if (geom.subset(math.index(i - 1, 1)) == Ly)
        {
            nf.subset(math.index(i - 1, 2), 0);
        }
    }
}
else if (end_index == 2)                // Fixed Ends
{
    for (var i = 1; i <= nnd; i++)
    {
        if (geom.subset(math.index(i - 1, 0)) == 0)
        {
            nf.subset(math.index(i - 1, 0), 0);
            nf.subset(math.index(i - 1, 1), 0);
            nf.subset(math.index(i - 1, 2), 0);
        }
        else if (geom.subset(math.index(i - 1, 1)) == 0)
        {
            nf.subset(math.index(i - 1, 0), 0);
            nf.subset(math.index(i - 1, 1), 0);
            nf.subset(math.index(i - 1, 2), 0);
        }
        else if (geom.subset(math.index(i - 1, 0)) == Lx)
        {
            nf.subset(math.index(i - 1, 1), 0);
        }
        else if (geom.subset(math.index(i - 1, 1)) == Ly)
        {
            nf.subset(math.index(i - 1, 2), 0);
        }
    }
}
else if(end_index == 3)                 // Beam supported slab
{
    // Write Code to implement the flexibility of beams.
}


// Count the free degrees of freedom

n = 0;

for (var i = 1; i <= nnd; i++)
{
    for (var j = 1; j <= nodof; j++)
    {
        if(nf.subset(math.index(i-1, j-1)) != 0)
        {
            n = n + 1;
            nf.subset(math.index(i-1, j-1), n);
        }
    }
}

//****************************************
// Load assignment


//Initialize the load matrix
load = math.zeros(nnd, 3);

var conc_load = 1000;           // This concentrated load will be updated from UI
var dead_load = density * thick;
var live_load = 5;              // To be updated from UI
var total_load_int = dead_load + live_load;         // Total intensity of loading
var total_load = total_load_int * Lx * Ly;          // Multiplied by area to get total load


if(loadtype_index == 1)                 // Concentrated Load
{
    for(var i = 1; i <= nnd; i++)
    {
        if (geom.subset(math.index(i - 1, 0)) == Lx && geom.subset(math.index(i - 1, 1)) == Ly) {
            load.subset(math.index(i - 1, 0), -1 * conc_load / 4);             // -1 to account for -Y axis and /4 as the conc load get divided among four symmetric plates
        }
    }
}
else if (loadtype_index == 2)           // UDL
{
    // Distribute to nodes based on the type of element
    if(etype_index == 1)
    {
        for(var i = 1; i <= nnd; i++)
        {
            count = 0;
            connec.forEach(function(value, index, matrix) {
                if(value == i) {
                    count = count + 1;
                }
            });
            load.subset(math.index(i - 1, 0), -1 * count * total_load / (4*nel));
        }
    }
    else if(etype_index == 2)
    {
        for(var i = 1; i <= nnd; i++)
        {
            count = 0;
            connec.forEach(function(value, index, matrix) {
                if(value == i) {
                    count = count + 1;
                }
            });
            load.subset(math.index(i - 1, 0), -1 * count * total_load / (8*nel));
        }
    }
}

// Define an inputs object that carries all the required inputs required for all the functions
var inputs = {nne: nne,
    nodof: nodof,
    eldof: eldof,
    geom: geom,
    connec: connec,
    nf: nf,
    dim: dim};

// ************************ End of Input *************

// Assemble the global force vector
var fg = math.zeros(n);

 for(var i = 1; i <= nnd; i++)
 {
    for (var j = 1; j <= nodof; j++)
    {
        if(nf.subset(math.index(i-1, j-1)) != 0)
        {
 //var row = nf.subset(math.index(i-1, j-1))
            fg.subset(math.index(nf.subset(math.index(i-1, j-1))-1), load.subset(math.index(i-1, j-1)));
        }
    }
 }

// Form the matrix containing the abscissa and weights of the Gauss points
sampb = gauss(ngpb);
samps = gauss(ngps);

// Numerical integration and assembly of global stiffness matrix
//
// Initialize the global stiffness matrix to zero
var kk = math.zeros(n,n);
var wi = 0;
var wj = 0;
var der = math.matrix();
var fun = math.matrix();
var jac = math.matrix();
var jacim = math.matrix();
var deriv = math.matrix();
var beeb = math.matrix();
var beebt = math.matrix();
var bees = math.matrix();
var beest = math.matrix();
var coord = math.matrix();
var g = math.matrix();

var d = 0;

// Initialize the element bending stiffness matrix to zero
var keb = math.zeros(eldof, eldof);

// Initialize the element shear stiffness matrix to zero
var kes = math.zeros(eldof, eldof);


for(var i = 1; i <= nel; i++)
{
    coord = platelem_q8(inputs, i).coord;
    g = platelem_q8(inputs, i).g;

    // Integrate element bending stiffness and assemble it in global matrix

    for(var ig = 1; ig <= ngpb; ig++)
    {
        wi = sampb.subset(math.index(ig-1, 1));
        for(var jg = 1; jg <= ngpb; jg++)
        {
            wj = sampb.subset(math.index(jg-1, 1));

            if(etype_index == 1)
            {
                der = fmlin(sampb, ig, jg).der;
                fun = fmlin(sampb, ig, jg).fun;
            }
            else if(etype_index == 2)
            {
                der = fmquad(sampb, ig, jg).der;
                fun = fmquad(sampb, ig, jg).fun;
            }

            jac = math.multiply(der, coord);
            d = math.det(jac);

            // Create a clone of the matrix, but a sylvester object now
            var jacs = Matrix.create(jac._data);

            // Calculate the inverse of the sylvester matrix
            var jacis = jacs.inverse();

            // Take it back to the mathjs matrix
            for(var p = 1; p <= jacis.rows(); p++)
            {
                for(var q = 1; q <= jacis.cols(); q++)
                {
                    jacim.subset(math.index(p-1,q-1), jacis.e(p,q)); // sylvester matrix brought to a mathjs matrix
                }
            }

            deriv = math.multiply(jacim, der);
            beeb = formbeeb(deriv, nne, eldof);

            // Tranpose the matrix
            var row_count = beeb._size[0];
            var col_count = beeb._size[1];
            for(var p = 0; p <= row_count-1; p++)
            {
                for(var q = 0; q <= col_count-1; q++)
                {
                    beebt.subset(math.index(q,p), beeb.subset(math.index(p, q)));
                }
            }

            // Integrate the stiffness matrix
            keb = math.add(keb, math.multiply(math.multiply(beebt, math.multiply(deeb, beeb)), d * wi * wj ));

        }
    }

    // Assemble global stiffness matrix
    kk = form_KK(inputs, kk, keb, g);

    // Integrate element Shear stiffness and assemble it in global matrix

    for(var ig = 1; ig <= ngps; ig++)
    {
        wi = sampb.subset(math.index(ig-1, 1));
        for(var jg = 1; jg <= ngps; jg++)
        {
            wj = sampb.subset(math.index(jg-1, 1));

            if(etype_index == 1)
            {
                der = fmlin(samps, ig, jg).der;
                fun = fmlin(samps, ig, jg).fun;
            }
            else if(etype_index == 2)
            {
                der = fmquad(samps, ig, jg).der;
                fun = fmquad(samps, ig, jg).fun;
            }

            jac = math.multiply(der, coord);
            d = math.det(jac);

            // Create a clone of the matrix, but a sylvester object now
            var jacs = Matrix.create(jac._data);

            // Calculate the inverse of the sylvester matrix
            var jacis = jacs.inverse();

            // Take it back to the mathjs matrix
            for(var p = 1; p <= jacis.rows(); p++)
            {
                for(var q = 1; q <= jacis.cols(); q++)
                {
                    jacim.subset(math.index(p-1,q-1), jacis.e(p,q)); // sylvester matrix brought to a mathjs matrix
                }
            }

            deriv = math.multiply(jacim, der);
            bees = formbees(deriv, fun, nne, eldof);

            // Tranpose the matrix
            var row_count = bees._size[0];
            var col_count = bees._size[1];
            for(var p = 0; p <= row_count-1; p++)
            {
                for(var q = 0; q <= col_count-1; q++)
                {
                    beest.subset(math.index(q,p), bees.subset(math.index(p,q)));
                }
            }

            // Integrate the stiffness matrix
            kes = math.add(kes, math.multiply(math.multiply(beest, math.multiply(dees, bees)),(5/6) * d * wi * wj));

        }
    }

    // Assemble global stiffness matrix
    kk = form_KK(inputs, kk, kes, g);
}

// **********************************************************
// ********* End of Assembly ********************************

var kks = Matrix.create(kk._data);
var fgs = Matrix.create(fg._data);
var delta = kks.solve(fgs);

console.log(nnd);
console.log(geom);
console.log(connec);
