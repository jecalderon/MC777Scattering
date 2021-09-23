// MC777Scattering.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

/********************************************
* Original work from Biomedical Optics Series
Steven L. Jacques, Scott A. Prahl 
 *  mc321.c    , in ANSI Standard C programing language
 *
 *  Monte Carlo simulation yielding spherical, cylindrical, and planar
 *    responses to an isotropic point source in an infinite homogeneous
 *    medium with no boundaries. This program is a minimal Monte Carlo
 *    program scoring photon distributions in spherical, cylindrical,
 *    and planar shells.
 *
 *  by Steven L. Jacques based on prior collaborative work
 *    with Lihong Wang, Scott Prahl, and Marleen Keijzer.
 *    partially funded by the NIH (R29-HL45045, 1991-1997) and
 *    the DOE (DE-FG05-91ER617226, DE-FG03-95ER61971, 1991-1999).
 *
 *  A published report illustrates use of the program:
 *    S. L. Jacques: "Light distributions from point, line, and plane
 *    sources for photochemical reactions and fluorescence in turbid
 *    biological tissues," Photochem. Photobiol. 67:23-32, 1998.
 *
 *  Trivial fixes to remove warnings SAP, 11/2017
 *  Theodore Info 07/2019: 
 *  The problem with the code was the fopen function. This function is deprecated
 *  and has been replace with the fopen_s function whose parameters are as follows
 *  fopen_s(<pointer to a file stream e.g FILE* >, <filename>, <options e.g w, r>)
 *  Modified by Abohweyere Oghenefejiro Theodore of Durham College Canada
 *  Modified by Jose E. Calderon University of Puerto Rico for a solution in VS C++ 2017
 *  
 *  There is a problem with the logic. If Absorption coeficient is smaller that 
 *  scattering coeeficient kills the simulation    8/6/2019
 *
 **********/

#include <math.h>
#include <stdio.h>
#include "pch.h"
#include <time.h>
#include <iostream>
#define Nbins	500
#define Nbinsp1	501
#define	PI          3.1415926
#define	LIGHTSPEED	2.997925E10 /* in vacuo speed of light [cm/s] */
#define ALIVE       1   		/* if photon not yet terminated */
#define DEAD        0    		/* if photon is to be terminated */
#define THRESHOLD   0.01		/* used in roulette */
#define CHANCE      0.1  		/* used in roulette */
#define COS90D      1.0E-6
 /* If cos(theta) <= COS90D, theta >= PI/2 - 1e-6 rad. */
#define ONE_MINUS_COSZERO 1.0E-12
	 /* If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */
	 /* If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */
#define SIGN(x)           ((x)>=0 ? 1:-1)
#define InitRandomGen    (double) RandomGen(0, 1, NULL)
	 /* Initializes the seed for the random number generator. */
#define RandomNum        (double) RandomGen(0, 1, NULL)
	 /* Calls for a random number from the randum number generator. */

/* DECLARE FUNCTION */
double RandomGen(char Type, long Seed, long *Status);
/* Random number generator */

inline void SphericalRadial(double x, double y, double z, double &r, double dr, short NR, short &ir){
	r = sqrt(x*x + y * y + z * z);    /* current spherical radial position */
	ir = (short)(r / dr);           /* ir = index to spatial bin */
	if (ir >= NR) ir = NR;        /* last bin is for overflow */
}

inline void CylindricalRadial(double x, double y, double &r, double dr, short NR, short &ir){
	r = sqrt(x*x + y * y);          /* current cylindrical radial position */
	ir = (short)(r / dr);           /* ir = index to spatial bin */
	if (ir >= NR) ir = NR;        /* last bin is for overflow */
}

inline void PlanarRadial(double z, double &r, double dr, short NR, short &ir){
	r = fabs(z);                  /* current planar radial position */
	ir = (short)(r / dr);           /* ir = index to spatial bin */
	if (ir >= NR) ir = NR;        /* last bin is for overflow */
}

int main() {

	/* Propagation parameters */
	register double		x, y, z;    /* photon position */
	register double		ux, uy, uz; /* photon trajectory as cosines */
	register double		uxx, uyy, uzz;	/* temporary values used during SPIN */
	register double		s;          /* step sizes. s = -log(RND)/mus [cm] */
	register double		costheta;   /* cos(theta) */
	register double		sintheta;   /* sin(theta) */
	register double		cospsi;     /* cos(psi) */
	register double		sinpsi;     /* sin(psi) */
	register double		psi;        /* azimuthal angle */
	register double		i_photon;   /* current photon */
	register double		W;          /* photon weight */
	register double		absorb;     /* weighted deposited in a step due to absorption */
	register short		photon_status;  /* flag = ALIVE=1 or DEAD=0 */

	/* other variables */
	register double	Csph[Nbinsp1];  /* spherical   photon concentration CC[ir=0..100] */
	register double	Ccyl[Nbinsp1];  /* cylindrical photon concentration CC[ir=0..100] */
	register double	Cpla[Nbinsp1];  /* planar      photon concentration CC[ir=0..100] */
	register double	Cobl[Nbinsp1];  /* Oblate      photon concentration CC[ir=0..100] */
	register double	Fsph;       /* fluence in spherical shell */
	register double	Fcyl;       /* fluence in cylindrical shell */
	register double	Fpla;       /* fluence in planar shell */
	register double	Fobl;       /* fluence in spheroidal shell */
	register double	mua;        /* absorption coefficient [cm^-1] */
	register double	mus;        /* scattering coefficient [cm^-1] */
	register double	g;          /* anisotropy [-] */
	register double	albedo;     /* albedo of tissue */
	register double	nt;         /* tissue index of refraction */
	register double	Nphotons;   /* number of photons in simulation */
	register short	NR;         /* number of radial positions */
	register double	radial_size;  /* maximum radial size */
	register double	r;          /* radial position */
	register double dr;         /* radial bin size */
	register short	ir;         /* index to radial position */
	register double  shellvolume;  /* volume of shell at radial position r */

	/* dummy variables */
	register double  rnd;        /* assigned random value 0-1 */
	register double	temp;    /* dummy variables */
	FILE*	target;     /* point to output file */

	clock_t tStart = clock();    /*  testing a time function  */

	/**** INPUT
	   Input the optical properties
	   Input the bin and array sizes
	   Input the number of photons
	*****/

	mua = 1.673;     /* cm^-1 */
	mus = 312.0;  /* cm^-1 */
	g = 0.9000;    /* The origina nummber is 0.9 */
	nt = 1.33;
	Nphotons = 10000; /* set number of photons in simulation */
	radial_size = 2.0;  /* cm, total range over which bins extend */
	NR = Nbins;	 /* set number of bins.  */
	   /* IF NR IS ALTERED, THEN USER MUST ALSO ALTER THE ARRAY DECLARATION TO A SIZE = NR + 1. */
	dr = radial_size / NR;  /* cm */
	albedo = mus / (mus + mua);


	/**** INITIALIZATIONS
	*****/
	i_photon = 0;
	InitRandomGen;
	for (ir = 0; ir <= NR; ir++) {
		Csph[ir] = 0;
		Ccyl[ir] = 0;
		Cpla[ir] = 0;
		Cobl[ir] = 0;  /*  Jose Calderon coordinate Oblate Spheroidal Coordinate */
	}

	/**** RUN
	   Launch N photons, initializing each one before progation.
	*****/
	do {


		/**** LAUNCH
		   Initialize photon position and trajectory.
		   Implements an isotropic point source.
		*****/
		i_photon += 1;	/* increment photon count */
		W = 1.0;                    /* set photon weight to one */
		photon_status = ALIVE;      /* Launch an ALIVE photon */

		x = 0;                      /* Set photon position to origin. */
		y = 0;
		z = 0;

		/* Randomly set photon trajectory to yield an isotropic source. */
		costheta = 2.0*RandomNum - 1.0;
		sintheta = sqrt(1.0 - costheta * costheta);	/* sintheta is always positive */
		psi = 2.0*PI*RandomNum;
		ux = sintheta * cos(psi);
		uy = sintheta * sin(psi);
		uz = costheta;


		/* HOP_DROP_SPIN_CHECK
		   Propagate one photon until it dies as determined by ROULETTE.
		*******/
		do {


			/**** HOP
			   Take step to new position
			   s = stepsize
			   ux, uy, uz are cosines of current photon trajectory
			*****/
			while ((rnd = RandomNum) <= 0.0);   /* yields 0 < rnd <= 1 */
			s = -log(rnd) / (mua + mus);          /* Step size.  Note: log() is base e */
			x += s * ux;                        /* Update positions. */
			y += s * uy;
			z += s * uz;


			/**** DROP
			   Drop photon weight (W) into local bin.
			*****/
			absorb = W * (1 - albedo);      /* photon weight absorbed at this step */
			W -= absorb;                  /* decrement WEIGHT by amount absorbed */

			/* spherical */
			SphericalRadial(x, y, z, r, dr, NR, ir);
			Csph[ir] += absorb;           /* DROP absorbed weight into bin */

			/* printf("Time taken Spheroida: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);     This time functio test*/

			/* cylindrical */
			CylindricalRadial(x, y, r, dr, NR, ir);
			Ccyl[ir] += absorb;           /* DROP absorbed weight into bin */

			/* planar */
			PlanarRadial(z, r, dr, NR, ir);
			Cpla[ir] += absorb;           /* DROP absorbed weight into bin */

			/* Oblate Spheroidal  */
			r = sqrt(x*x / (1 + sqrt(2.0)) - y * y + z * z);  /* current spheroidal radial position wher focal point equal minor axis */
			ir = (short)(r / dr);           /* ir = index to spatial bin */
			if (ir >= NR) ir = NR;			/* last bin is for overflow */
			Cobl[ir] += absorb;				/* DROP absorbed weight into bin */


		 /**** SPIN
			Scatter photon into new trajectory defined by theta and psi.
			Theta is specified by cos(theta), which is determined
			based on the Henyey-Greenstein scattering function.
			Convert theta and psi into cosines ux, uy, uz.
		 *****/
		 /* Sample for costheta */
			rnd = RandomNum;
			if (g == 0.0)
				costheta = 2.0*rnd - 1.0;
			else {
				double temp = (1.0 - g * g) / (1.0 - g + 2 * g*rnd);
				costheta = (1.0 + g * g - temp * temp) / (2.0*g);
			}
			sintheta = sqrt(1.0 - costheta * costheta); /* sqrt() is faster than sin(). */

			/* Sample psi. */
			psi = 2.0*PI*RandomNum;
			cospsi = cos(psi);
			if (psi < PI)
				sinpsi = sqrt(1.0 - cospsi * cospsi);     /* sqrt() is faster than sin(). */
			else
				sinpsi = -sqrt(1.0 - cospsi * cospsi);

			/* New trajectory. */
			if (1 - fabs(uz) <= ONE_MINUS_COSZERO) {      /* close to perpendicular. */
				uxx = sintheta * cospsi;
				uyy = sintheta * sinpsi;
				uzz = costheta * SIGN(uz);   /* SIGN() is faster than division. */
			}
			else {					/* usually use this option */
				temp = sqrt(1.0 - uz * uz);
				uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta;
				uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta;
				uzz = -sintheta * cospsi * temp + uz * costheta;
			}

			/* Update trajectory */
			ux = uxx;
			uy = uyy;
			uz = uzz;


			/**** CHECK ROULETTE
			   If photon weight below THRESHOLD, then terminate photon using Roulette technique.
			   Photon has CHANCE probability of having its weight increased by factor of 1/CHANCE,
			   and 1-CHANCE probability of terminating.
			*****/
			if (W < THRESHOLD) {
				if (RandomNum <= CHANCE)
					W /= CHANCE;
				else photon_status = DEAD;
			}


		} /* end STEP_CHECK_HOP_SPIN */
		while (photon_status == ALIVE);
		/* If photon dead, then launch new photon. */
	} /* end RUN */
	while (i_photon < Nphotons);


	/**** SAVE
	   Convert data to relative fluence rate [cm^-2] and save to file called "mcmin321.out".
	   Theodore: Here fopen(<filename>, <options>) has been replaced with
	   fopen_s(<pointer to a file stream e.g FILE* >, <filename>, <options e.g w, r>)
	*****/
    fopen_s(&target, "mc321_.out", "w");

	/* print header */
	fprintf(target, "number of photons = %f\n", Nphotons);
	fprintf(target, "bin size = %5.5f [cm] \n", dr);
	fprintf(target, "last row is overflow. Ignore.\n");

	/* print column titles */
	fprintf(target, "r [cm] \t Fsph [1/cm2] \t Fcyl [1/cm2] \t Fpla [1/cm2] \t Fobl [1/cm2]\n");

	/* print data:  radial position, fluence rates for 3D, 2D, 1D geometries */
	for (ir = 0; ir <= NR; ir++) {
		/* r = sqrt(1.0/3 - (ir+1) + (ir+1)*(ir+1))*dr; */
		r = (ir + 0.5)*dr;
		shellvolume = 4.0*PI*r*r*dr; /* per spherical shell */
		Fsph = Csph[ir] / Nphotons / shellvolume / mua;
		shellvolume = 2.0*PI*r*dr;   /* per cm length of cylinder */
		Fcyl = Ccyl[ir] / Nphotons / shellvolume / mua;
		shellvolume = dr;            /* per cm2 area of plane */
		Fpla = Cpla[ir] / Nphotons / shellvolume / mua;
		shellvolume = 2.0*(1 + sqrt(2.0))*PI*r*r*dr; /* per spheroidal shell */
		Fobl = Cobl[ir] / Nphotons / shellvolume / mua;
		fprintf(target, "%5.5f \t %4.3e \t %4.3e \t %4.3e \t %4.3e \n", r, Fsph, Fcyl, Fpla, Fobl);
	}

	fclose(target);
	
	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	

} /* end of main */



/* SUBROUTINES */

/**************************************************************************
 *	RandomGen
 *      A random number generator that generates uniformly
 *      distributed random numbers between 0 and 1 inclusive.
 *      The algorithm is based on:
 *      W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.
 *      Flannery, "Numerical Recipes in C," Cambridge University
 *      Press, 2nd edition, (1992).
 *      and
 *      D.E. Knuth, "Seminumerical Algorithms," 2nd edition, vol. 2
 *      of "The Art of Computer Programming", Addison-Wesley, (1981).
 *
 *      When Type is 0, sets Seed as the seed. Make sure 0<Seed<32000.
 *      When Type is 1, returns a random number.
 *      When Type is 2, gets the status of the generator.
 *      When Type is 3, restores the status of the generator.
 *
 *      The status of the generator is represented by Status[0..56].
 *
 *      Make sure you initialize the seed before you get random
 *      numbers.
 ****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

double RandomGen(char Type, long Seed, long *Status) {
	static long i1, i2, ma[56];   /* ma[0] is not used. */
	long        mj, mk;
	short       i, ii;

	if (Type == 0) {              /* set seed. */
		mj = MSEED - (Seed < 0 ? -Seed : Seed);
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for (i = 1; i <= 54; i++) {
			ii = (21 * i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if (mk < MZ)
				mk += MBIG;
			mj = ma[ii];
		}
		for (ii = 1; ii <= 4; ii++)
			for (i = 1; i <= 55; i++) {
				ma[i] -= ma[1 + (i + 30) % 55];
				if (ma[i] < MZ)
					ma[i] += MBIG;
			}
		i1 = 0;
		i2 = 31;
	}
	else if (Type == 1) {       /* get a number. */
		if (++i1 == 56)
			i1 = 1;
		if (++i2 == 56)
			i2 = 1;
		mj = ma[i1] - ma[i2];
		if (mj < MZ)
			mj += MBIG;
		ma[i1] = mj;
		return (mj * FAC);
	}
	else if (Type == 2) {       /* get status. */
		for (i = 0; i < 55; i++)
			Status[i] = ma[i + 1];
		Status[55] = i1;
		Status[56] = i2;
	}
	else if (Type == 3) {       /* restore status. */
		for (i = 0; i < 55; i++)
			ma[i + 1] = Status[i];
		i1 = Status[55];
		i2 = Status[56];
	}
	else
		puts("Wrong parameter to RandomGen().");
	return (0);
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
