/* Generates the incident beam
 *
 * Lminus beam is based on: G. Gouesbet, B. Maheu, G. Grehan, "Light scattering from a sphere arbitrary located
 * in a Gaussian beam, using a Bromwhich formulation", J.Opt.Soc.Am.A 5,1427-1443 (1988).
 * Eq.(22) - complex conjugate.
 *
 * Davis beam is based on: L. W. Davis, "Theory of electromagnetic beams," Phys.Rev.A 19, 1177-1179 (1979).
 * Eqs.(15a),(15b) - complex conjugate; in (15a) "Q" changed to "Q^2" (typo).
 *
 * Barton beam is based on: J. P. Barton and D. R. Alexander, "Fifth-order corrected electromagnetic-field components
 * for a fundamental Gaussian-beam," J.Appl.Phys. 66,2800-2802 (1989).
 * Eqs.(25)-(28) - complex conjugate.
 *
 * Copyright (C) ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#include "const.h" // keep this first
// project headers
#include "cmplx.h"
#include "comm.h"
#include "io.h"
#include "interaction.h"
#include "param.h"
#include "vars.h"
// system headers
#include <stdio.h>
#include <stdlib.h> // for abs()
#include <string.h>

#include <math.h>

// SEMI-GLOBAL VARIABLES

// defined and initialized in param.c
extern const int beam_Npars;
extern const double beam_pars[];
extern const char *beam_fnameY;
extern const char *beam_fnameX;
extern const opt_index opt_beam;

// used in CalculateE.c
double C0dipole,C0dipole_refl; // inherent cross sections of exciting dipole (in free space and addition due to surface)
int vorticity;                 // Vorticity of vortex beams (besN for Bessel beams)
// used in param.c
const char *beam_descr; // string for log file with beam parameters
/* Propagation (phase) directions of secondary incident beams above (refl) and below (tran) the surface (unit vectors)
 * When msub is complex, one of this doesn't tell the complete story, since the corresponding wave is inhomogeneous,
 * given by the complex wavenumber ktVec
 */
double prIncRefl[3],prIncTran[3];

// LOCAL VARIABLES
static double s,s2;            // beam confinement factor and its square
static double scale_x,scale_z; // multipliers for scaling coordinates
static doublecomplex ki,kt;    // abs of normal components of k_inc/k0, and ktran/k0
static doublecomplex ktVec[3]; // k_tran/k0
static double p0;              // amplitude of the incident dipole moment
#ifndef NO_FORTRAN
static int besN;                  // Bessel beam order
static double besAlpha;           // half-cone angle (in radians)
static doublecomplex besKt,besKz; // wave-vector components (transverse and longitudinal)
static doublecomplex besM[4];     // Matrix M defining the generalized Bessel beam
#endif
static int l;		//Radial order for Laguerre-Gaussian beam
static int p;		//Azimuthal order for Laguerre-Gaussian beam
/* TO ADD NEW BEAM
 * Add here all internal variables (beam parameters), which you initialize in InitBeam() and use in GenerateB()
 * afterwards. If you need local, intermediate variables, put them into the beginning of the corresponding function.
 * Add descriptive comments, use 'static'.
 */

// EXTERNAL FUNCTIONS

#ifndef NO_FORTRAN
void bjndd_(const int *n,const double *x,double *bj,double *dj,double *fj);
#endif

//======================================================================================================================

void InitBeam(void)
// initialize beam; produce description string
{
	double w0; // beam width
	const char *tmp_str; // temporary string
	/* TO ADD NEW BEAM
	 * Add here all intermediate variables, which are used only inside this function.
	 */

	// initialization of global option index for error messages
	opt=opt_beam;
	// beam initialization
	beam_asym=(beam_center_0[0]!=0 || beam_center_0[1]!=0 || beam_center_0[2]!=0);
	if (!beam_asym) vInit(beam_center); // not to calculate it for each orientation
	vorticity = 0;
	switch (beamtype) {
		case B_PLANE:
			if (IFROOT) beam_descr="plane wave";
			if (surface) {
				/* here we assume that prop_0 will not change further (e.g., by rotation of particle),
				 * i.e. prop=prop_0 in GenerateBeam() below
				 */
				if (prop_0[2]==0) PrintError("Ambiguous setting of beam propagating along the surface. Please specify "
					"the incident direction to have (arbitrary) small positive or negative z-component");
				if (msubInf && prop_0[2]>0) PrintError("Perfectly reflecting surface ('-surf ... inf') is incompatible "
					"with incident direction from below (including the default one)");
				// Here we set ki,kt,ktVec and propagation directions prIncRefl,prIncTran
				if (prop_0[2]>0) { // beam comes from the substrate (below)
					// here msub should always be defined
					inc_scale=1/creal(msub);
					ki=msub*prop_0[2];
					/* Special case for msub near 1 to remove discontinuities for near-grazing incidence. The details
					 * are discussed in CalcFieldSurf() in crosssec.c.
					 */
					if (cabs(msub-1)<ROUND_ERR && cabs(ki)<SQRT_RND_ERR) kt=ki;
					else kt=cSqrtCut(1 - msub*msub*(prop_0[0]*prop_0[0]+prop_0[1]*prop_0[1]));
					// determine propagation direction and full wavevector of wave transmitted into substrate
					ktVec[0]=msub*prop_0[0];
					ktVec[1]=msub*prop_0[1];
					ktVec[2]=kt;
				}
				else if (prop_0[2]<0) { // beam comes from above the substrate
					inc_scale=1;
					ki=-prop_0[2]; // always real
					if (!msubInf) {
						// same special case as above
						if (cabs(msub-1)<ROUND_ERR && cabs(ki)<SQRT_RND_ERR) kt=ki;
						else kt=cSqrtCut(msub*msub - (prop_0[0]*prop_0[0]+prop_0[1]*prop_0[1]));
						// determine propagation direction of wave transmitted into substrate
						ktVec[0]=prop_0[0];
						ktVec[1]=prop_0[1];
						ktVec[2]=-kt;
					}
				}
				else {
					LogError(ONE_POS,"Ambiguous setting of beam propagating along the surface. Please specify the"
					"incident direction to have (arbitrary) small positive or negative z-component");
				}
				vRefl(prop_0,prIncRefl);
				if (!msubInf) {
					vReal(ktVec,prIncTran);
					vNormalize(prIncTran);
				}
			}
			printf("LOG PRINT; ---- 	OUTSIDE surf ---- \n");
			return;
		case B_DIPOLE:
			if (surface) {
				if (beam_center_0[2]<=-hsub)
					PrintErrorHelp("External dipole should be placed strictly above the surface");
				inc_scale=1; // but scaling of Mueller matrix is weird anyway
			}
			/* definition of p0 is important for scaling of many scattering quantities (that are normalized to incident
			 * irradiance). Alternative definition is p0=1, but then the results will scale with unit of length
			 * (breaking scale invariance)
			 */
			p0=1/(WaveNum*WaveNum*WaveNum);
			if (IFROOT) beam_descr="point dipole";
			return;
		case B_LAGUERRE:
			l=beam_pars[0];
			p=beam_pars[1];
			vorticity=p;
			//####################################
			//## POSSIBLY CHANGE VORTICITY HERE ##
			//####################################
			//if (IFROOT) {
			//	tmp_str="Laguerre-Gaussian formulation";
			//	break;
			//}
			return;
		case B_LMINUS:
		case B_DAVIS3:
		case B_BARTON5:
			if (surface) PrintError("Currently, Gaussian incident beam is not supported for '-surf'");
			// initialize parameters
			w0=beam_pars[0];
			TestPositive(w0,"beam width");
			double wavelength = 2.0*PI/(WaveNum);
			s=1/(WaveNum*w0);
			s2=s*s;
			scale_x=1/w0;
			scale_z=s*scale_x; // 1/(k*w0^2)
			// beam info
			if (IFROOT) {
				switch (beamtype) {
					case B_LMINUS:
						tmp_str="L- approximation";
						break;
					case B_DAVIS3:
						tmp_str="3rd order approximation, by Davis";
						break;
					case B_BARTON5:
						tmp_str="5th order approximation, by Barton";
						break;
					default: LogError(ONE_POS,"Incompatibility error in GenerateB");
				}
				beam_descr=dyn_sprintf("Gaussian beam (%s)\n"
				                       "\tWidth="GFORMDEF" (confinement factor s="GFORMDEF")",tmp_str,w0,s);
			}
			return;
#ifndef NO_FORTRAN
		case B_BES_CS:
		case B_BES_CSp:
		case B_BES_M:
		case B_BES_LE:
		case B_BES_LM:
		case B_BES_TEL:
		case B_BES_TML:
			if (surface) PrintError("Currently, Bessel incident beam is not supported for '-surf'");
			// initialize parameters
			ConvertToInteger(beam_pars[0],"beam order",&besN);
			TestRangeII(besN,"beam order (might cause the incorrect calculation of Bessel function)",-50,50);
			vorticity=besN;
			besAlpha=Deg2Rad(beam_pars[1]);
			besKt=WaveNum*sin(besAlpha);
			besKz=WaveNum*cos(besAlpha);
			/* redundant initialization to remove warnings. Compiler is not sure about constancy of IFROOT due to
			 * ringid being a global variable
			 */
			tmp_str="";
			switch (beamtype) { // definition of elements of matrix M ((Mex,Mey),(Mmx,Mmy))
				case B_BES_CS:
					printf("LOG PRINT; ---- Inside Bessel CS Beam SPECIFIC ---- \n");
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					besM[0]=0.5;  besM[1]=0;
					besM[2]=0;    besM[3]=0.5;
					if (IFROOT) tmp_str="circularly symmetric energy density";
					break;
				case B_BES_CSp:
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					besM[0]=0.5;  besM[1]=0;
					besM[2]=0;    besM[3]=-0.5;
					if (IFROOT) tmp_str="circularly symmetric energy density, alternative type";
					break;
				case B_BES_M:
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					if (beam_Npars==6) {
						besM[0]=beam_pars[2];
						besM[1]=beam_pars[3];
						besM[2]=beam_pars[4];
						besM[3]=beam_pars[5];
					}
					else {
						besM[0]=beam_pars[2]+I*beam_pars[6];
						besM[1]=beam_pars[3]+I*beam_pars[7];
						besM[2]=beam_pars[4]+I*beam_pars[8];
						besM[3]=beam_pars[5]+I*beam_pars[9];
					}
					if (IFROOT) tmp_str="generalized";
					break;
				case B_BES_LE:
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					besM[0]=0;  besM[1]=0;
					besM[2]=0;  besM[3]=1;
					if (IFROOT) tmp_str="linear electric field";
					break;
				case B_BES_LM:
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					besM[0]=0;  besM[1]=1;
					besM[2]=0;  besM[3]=0;
					if (IFROOT) tmp_str="linear magnetic field";
					break;
				/* TODO: for the following two types, both 0 and 90 degrees should be fine, but may require some
				 * rearrangement of formulae. Then the tests for range of alpha should be moved outside of the case
				 */
				case B_BES_TEL:
					TestRangeNN(beam_pars[1],"half-cone angle for TEL type",0,90);
					besM[0]=-WaveNum/besKt;  besM[1]=0;
					besM[2]=0;               besM[3]=besKz/besKt;
					if (IFROOT) tmp_str="linear component of the TE";
					break;
				case B_BES_TML:
					TestRangeNN(beam_pars[1],"half-cone angle for TML type",0,90);
					besM[0]=0;              besM[1]=besKz/besKt;
					besM[2]=WaveNum/besKt;  besM[3]=0;
					if (IFROOT) tmp_str="linear component of the TM";
					break;
				default: LogError(ONE_POS,"Incompatibility error in GenerateB");
			}
			// TODO: some symmetries can be retained in some special cases
			symR=symX=symY=symZ=false;
			// beam info
			if (IFROOT) beam_descr=dyn_sprintf("Bessel beam (%s)\n"
				                               "\tOrder: %d, half-cone angle: "GFORMDEF" deg",
				                               tmp_str,besN,beam_pars[1]);
			return;
#endif // !NO_FORTRAN
		case B_READ:
			// the safest is to assume cancellation of all symmetries
			symX=symY=symZ=symR=false;
			if (surface) inc_scale=1; // since we can't know it, we assume the default case
			if (IFROOT) {
				if (beam_Npars==1) beam_descr=dyn_sprintf("specified by file '%s'",beam_fnameY);
				else beam_descr=dyn_sprintf("specified by files '%s' and '%s'",beam_fnameY,beam_fnameX);
			}
			// this type is unaffected by beam_center
			return;
	}
	LogError(ONE_POS,"Unknown type of incident beam (%d)",(int)beamtype);
	/* TO ADD NEW BEAM
	 * add a case above. Identifier ('B_...') should be defined inside 'enum beam' in const.h. The case should
	 * 1) save all the input parameters from array 'beam_pars' to local variables (defined in the beginning of this
	 *    source files)
	 * 2) test all input parameters (for that you're encouraged to use functions from param.h since they would
	 *    automatically produce informative output in case of error).
	 * 3) the symmetry breaking due to prop or beam_center is taken care of in VariablesInterconnect() in param.c.
	 *    But if there are other reasons why beam would break any symmetry, corresponding variable should be set to
	 *    false here. Do not set any of them to true, as they can be set to false by other factors.
	 *    symX, symY, symZ - symmetries of reflection over planes YZ, XZ, XY respectively.
	 *    symR - symmetry of rotation for 90 degrees over the Z axis
	 * 4) initialize the following:
	 *    beam_descr - descriptive string, which will appear in log file  (should NOT end with \n).
	 *    vorticity - (only for vortex beam) integer value, how many turns the phase experience, when one makes a full
	 *                turn around the beam axis.
	 * 5) Consider the case of surface (substrate near the particle). If the new beam type is incompatible with it, add
	 *    an explicit exception, like "if (surface) PrintErrorHelp(...);". Otherwise, you also need to define inc_scale.
	 * All other auxiliary variables, which are used in beam generation (GenerateB(), see below), should be defined in
	 * the beginning of this file. If you need temporary local variables (which are used only in this part of the code),
	 * define them in the beginning of this function.
	 */
}


//======================================================================================================================


doublecomplex integrate_simpsons_rule(double limit_lower, double limit_upper, int samples, doublecomplex (*cFunction)(double, int, double, double, double, double, int, int, double, double, doublecomplex, doublecomplex), int component, double rho, double phi, double z, double k, int l, int p, double z_R, double w0, doublecomplex alpha, doublecomplex beta) {
	/*
	. Evaluates an integral using simpson's rule
	*/
	int i_ind;
	doublecomplex integral = 0.0 +0.0*I;
	double sample_step = (limit_upper-limit_lower)/(samples);

	//Goes to 0 at end points (problems if included)
	//integral += cFunction( (limit_lower), component, rho, phi, z, k, l, p, z_R, w0, alpha, beta);
	for(int i_ind=1; i_ind<samples; i_ind++) {
		//Loop and evaluate at the same time
		if(i_ind%2==0) {
			//2* for evens
			integral += 2.0*cFunction( (limit_lower +i_ind*sample_step), component, rho, phi, z, k, l, p, z_R, w0, alpha, beta);
		} else {
			//4* for odds
			integral += 4.0*cFunction( (limit_lower +i_ind*sample_step), component, rho, phi, z, k, l, p, z_R, w0, alpha, beta);
		}
	}
	//Goes to 0 at end points (problems if included)
	//integral += cFunction( (limit_upper), component, rho, phi, z, k, l, p, z_R, w0, alpha, beta);
	return integral*(sample_step/3.0);
}
double laguerre_getParam_w(double zeta, double w0) {
	/*
	. Get the 'w' term used in laguerre-gaussian beam
	*/
	return w0*sqrt(1.0+(zeta*zeta));
}
doublecomplex laguerre_getParam_Llp(doublecomplex value, int l, int p) {
	/*
	. Get associated legendre polynomial of given order
	. This matches the values given by python clpmn() in the SPECIFIC l=0 case => associated 
		legendre polynomials are equivalent to legendre polynomials

	****************************************************************************************
	. TEMPORARY FIX; Currently the associated legendre polynomial produced ONLY corresponds
		to values l=0, p=8
	****************************************************************************************
	*/
	doublecomplex L = 0.0 +0.0*I;	//Generalised implementation ### UNFINISHED PLACEHOLDER ###
	doublecomplex L_0_8 = ( (6435.0)*cpow(value,8) -(12012.0)*cpow(value,6) +(6930.0)*cpow(value,4) -(1260.0)*cpow(value,2) +(35.0) )/128.0;
	return L_0_8;
}
doublecomplex laguerre_getParam_Bessel(doublecomplex value, int l) {
	/*
	. Get bessel function
	#####################################################
	## MAY ALSO NOT BE CORRECT FORM FOR COMPLEX VALUES ##
	#####################################################
	*/
	int terms = 6;	//# of terms to take in Taylor series
	doublecomplex bessel_value = 0.0 +0.0*I;
	for(int i_ind=0; i_ind<terms; i_ind++) {
		doublecomplex set_1 = (pow(-1,i_ind))/(factorial(i_ind)*factorial(i_ind+l));
		doublecomplex set_2 = cpow( (value/2.0), (2*i_ind +l));
		bessel_value += set_1*set_2;
	}

	//###############################
	//## VALUE SEEMS EXTREMELY LOW ##
	//###############################
	//bessel_value = 1.0 +0.0*I;
	//printf("bessel= %e, %e \n", creal(bessel_value), cimag(bessel_value));
	return bessel_value;
}
doublecomplex laguerre_getParam_C(int l, int p) {
	/*
	. Get normalisation constant C, used in u_l_p function
	. From “Gaussian Beams In Optics of Course” paper
	*/
	//#####
	//## DIFFERENT FORMULATIONS MAY HAVE COMPLEX PARTS --> HENCE BETTER TO ASSUME COMPLEX
	//#####
	//(factorial(p))*sqrt( (2)/(PI*factorial(p)*( factorial( abs(l) + p ) )) );
	doublecomplex C_formualtion = (factorial(p))*sqrt( (2)/(PI*factorial(p)*( factorial( abs(l) + p ) )) );
	return C_formualtion;
}
doublecomplex laguerre_getParam_Ulp(doublecomplex C, int l, int p, double rho, double phi, double k, double z, double z_R, double w0) {
	/*
	. Get u_l_p, which describes behaviour of laguerre-gaussian beam
	*/
	//#####
	//## SOME INLINE COMPLEX VALUES HERE --> THINK IS FINE BUT WORTH CHECKING
	//#####
	doublecomplex set_1 = (C)/csqrt(1 +((z*z)/(z_R*z_R)));
	doublecomplex set_2 = cpow( (rho*csqrt(2))/(laguerre_getParam_w(z, w0)), l);
	doublecomplex set_3 = laguerre_getParam_Llp((2.0*(rho*rho))/(cpow(laguerre_getParam_w(z, w0),2)), l, p);
	doublecomplex set_4 = imExp( (-(rho*rho))/(cpow(laguerre_getParam_w(z, w0),2)) );	//<------ This term gets very small, ends up setting E to 0
	doublecomplex set_5 = imExp( (-I*k*(rho*rho)*z)/(2.0*((z*z)+(z_R*z_R))) );
	doublecomplex set_6 = imExp(I*l*phi);
	doublecomplex set_7 = imExp( (I*(2.0*p +l+1))*(atan( (z)/(z_R) )) );
	return set_1*set_2*set_3*set_4*set_5*set_6*set_7;
}
doublecomplex laguerre_getParam_E(double kappa, double rho, double phi, double z, double k, int l, int p, double z_R, double w0) {
	/*
	. Get E magnitude term for laguerre-gaussian beam
	*/
	doublecomplex C = laguerre_getParam_C(l, p);//1.0 +0.0*I;
	doublecomplex set_1 = laguerre_getParam_Ulp(C, l, p, rho, phi, k, z, z_R, w0);
	doublecomplex set_2 = imExp( (-k*(kappa*kappa)*z_R)/(2.0*((k*k) - (kappa*kappa))) );
	doublecomplex set_3 = cpow( ((kappa*kappa))/((k*k) - (kappa*kappa)), (2.0*p+l+1.0)/(2.0) );
	doublecomplex set_4 = csqrt( ((k*k))/((k*k) - (kappa*kappa)) );
	doublecomplex set_5 = set_1*set_2*set_3*set_4;
	return set_5;
}
doublecomplex laguerre_getParam_E_component(double kappa, int component, double rho, double phi, double z, double k, int l, int p, double z_R, double w0, doublecomplex alpha, doublecomplex beta) {
	/*
	. Returns magnitude of the ith component of E (x,y or z)
	*/
	doublecomplex E = laguerre_getParam_E(kappa, rho, phi, z, k, l, p, z_R, w0);
	if(component == 0) {
		doublecomplex set_1 = imExp(I*l*phi);
		doublecomplex set_2 = imExp(I*z*csqrt((k*k)-(kappa*kappa)));
		doublecomplex set_3 = ( alpha*laguerre_getParam_Bessel(kappa*rho, l) );
		//doublecomplex set_3 = ( alpha*(1.0) );
		//printf("E = %e, %e \n",creal(E), cimag(E));
		//printf("set_1 = %e, %e \n",creal(set_1), cimag(set_1));
		//printf("set_2 = %e, %e \n",creal(set_2), cimag(set_2));
		//printf("set_3 = %e, %e \n",creal(set_3), cimag(set_3));
		return E*set_1*set_2*set_3;
	} else if(component == 1) {
		//Same as X component, but beta NOT alpha term
		double temp_v0 = creal(kappa)*rho;
		double temp_v1, temp_v2, temp_v3;
		//bjndd_(&l, &temp_v0, &temp_v1, &temp_v2, &temp_v3 );		//Bessel function, only care about t2, not t3 or 4	//##### CAREFUL OF TYPES FOR THIS RETURN #####
		doublecomplex set_1 = imExp(I*l*phi);
		doublecomplex set_2 = imExp(I*z*csqrt((k*k)-(kappa*kappa)));
		doublecomplex set_3 = ( beta*laguerre_getParam_Bessel(kappa*rho, l) );
		return E*set_1*set_2*set_3;
	} else if(component == 2) {
		doublecomplex set_1 = imExp(I*l*phi);
		doublecomplex set_2 = imExp(I*z*csqrt((k*k)-(kappa*kappa)));
		doublecomplex set_3 = (kappa)/(2.0*csqrt((k*k) - (kappa*kappa)));
		doublecomplex set_4 = (I*alpha - beta)*(imExp(-I*phi))*laguerre_getParam_Bessel(kappa*rho, l-1) - (I*alpha + beta)*(imExp(I*phi))*laguerre_getParam_Bessel(kappa*rho, l+1);
		return E*set_1*set_2*set_3*set_4;
	} else {
		printf("Unrecognised component \n");
		return 0.0 +0.0*I;
	}
}
int factorial(int value) {
	/*
	. Returns the factorial of a positive integer
	*/
	int break_condition = 50;	//Break after 50 loops to prevent inf. loop, factorials deal with here never more than ~10
	int factorial_value = 1;
	if(value > 0) {
		//Returns 1 if 0 or less
		while(value > 0) {
			factorial_value *= value;
			value = value-1;
			break_condition -= 1;
			if(break_condition < 0) {
				break;}
		}
	}
	return factorial_value;
}


//======================================================================================================================


void GenerateB (const enum incpol which,   // x - or y polarized incident light
                doublecomplex *restrict b) // the b vector for the incident field
// generates incident beam at every dipole
{
	size_t i,j;
	doublecomplex psi0,Q,Q2;
	doublecomplex v1[3],v2[3],v3[3],gt[6];
	double ro2,ro4;
	double x,y,z,x2_s,xy_s;
	doublecomplex t1,t2,t3,t4,t5,t6,t7,t8,ctemp;	//Temporary doublecomplex
	const double *ex; // coordinate axis of the beam reference frame
	double ey[3];
	double r1[3];
	/* complex wave amplitudes of transmitted wave (with phase relative to beam center);
	 * The transmitted wave can be inhomogeneous wave (when msub is complex), then eIncTran (e) is normalized
	 * counter-intuitively. Before multiplying by tc, it satisfies (e,e)=1!=||e||^2. This normalization is consistent
	 * with used formulae for transmission coefficients. So this transmission coefficient is not (generally) equal to
	 * the ratio of amplitudes of the electric fields. In particular, when E=E0*e, ||E||!=|E0|*||e||, where
	 * ||e||^2=(e,e*)=|e_x|^2+|e_y|^2+|e_z|^2=1
	 */
	doublecomplex eIncTran[3];
#ifndef NO_FORTRAN
	// for Bessel beams
	int n1,q;
	doublecomplex vort;  // vortex phase of Bessel beam rotated by 90 deg
	doublecomplex fn[5]; // for general functions f(n,ro,phi) of Bessel beams (fn-2, fn-1, fn, fn+1, fn+2, respectively)
	double phi,arg,td1[abs(besN)+3],td2[abs(besN)+3],jn1[abs(besN)+3]; // for Bessel beams
#endif
	const char *fname;

	/* TO ADD NEW BEAM
	 * Add here all intermediate variables, which are used only inside this function. You may as well use 't1'-'t8'
	 * variables defined above.
	 */

	// set reference frame of the beam; ez=prop, ex - incident polarization
	if (which==INCPOL_Y) {
		ex=incPolY;
		vMultScal(-1,incPolX,ey);
	}
	else { // which==INCPOL_X
		ex=incPolX;
		vCopy(incPolY,ey);
	}

	switch (beamtype) {
		case B_PLANE: // plane is separate to be fast (for non-surface)
			if (surface) {
				/* With respect to normalization we use here the same assumption as in the free space - the origin is in
				 * the particle center, and amplitude of incoming plane wave is equal to 1. Then irradiance of the beam
				 * coming from below is c*Re(msub)/(8pi), different from that coming from above.
				 * Original incident (incoming) beam propagating from the vacuum (above) is Exp(i*k*r.a), while - from
				 * the substrate (below) is Exp(i*k*msub*r.a). We assume that the incoming beam is homogeneous in its
				 * original medium.
				 */
				double hbeam=hsub+beam_center[2]; // height of beam center above the surface
				if (prop[2]>0) { // beam comes from the substrate (below)
					doublecomplex tc; // transmission coefficients
					//  determine amplitude of the transmitted wave; here msub is always defined
					if (which==INCPOL_Y) { // s-polarized
						cvBuildRe(ex,eIncTran);
						tc=FresnelTS(ki,kt);
					}
					else { // p-polarized
						crCrossProd(ey,ktVec,eIncTran);
						tc=FresnelTP(ki,kt,1/msub);
					}
					// phase shift due to the beam center relative to the origin and surface
					cvMultScal_cmplx(tc*cexp(I*WaveNum*((kt-ki)*hbeam-crDotProd(ktVec,beam_center))),eIncTran,eIncTran);
					// main part
					for (i=0;i<local_nvoid_Ndip;i++) {
						j=3*i;
						// b[i] = eIncTran*exp(ik*kt.r); subtraction of beam_center is avoided by the factor above
						cvMultScal_cmplx(cexp(I*WaveNum*crDotProd(ktVec,DipoleCoord+j)),eIncTran,b+j);
					}
				}
				else if (prop[2]<0) { // beam comes from above the substrate
					/* The following code takes extra care to be stable (not to lose precision) for grazing incidence.
					 * While it seems more complicated for general incidence, it requires only a few more complex
					 * arithmetic operations that are negligible compared to two complex exponents. For grazing
					 * incidence, it is not only stable but may also be faster than the straightforward computation,
					 * since one of the complex exponents is computed from a very small argument
					 */
					// determine reflection coefficient + 1
					doublecomplex rcp1;
					if (which==INCPOL_Y) { // s-polarized
						if (msubInf) rcp1=0;
						else rcp1=FresnelTS(ki,kt);
					}
					else { // p-polarized
						if (msubInf) rcp1=2;
						else rcp1=msub*FresnelTP(ki,kt,msub);
					}
					// main part
					for (i=0;i<local_nvoid_Ndip;i++) {
						j=3*i;
						vSubtr(DipoleCoord+j,beam_center,r1); // beam_center affects only the common phase factor
						/* b[i] = ex*exp(ik*r.a) + eIncRefl*exp[ik(prIncRefl.r-2az*hbeam)],
						 * where prIncRefl is prop with reflected z-component, while eIncRefl=rc*ex for s-polarization
						 * and (additionally) with reflected transverse (x,y) components for p-polarization
						 */
						ctemp=imExp(WaveNum*DotProd(r1,prop)); // exp(ik*r.a)
						t1=imExpM1(-2*WaveNum*prop[2]*(r1[2]+hbeam));
						/* t2=exp(ik*r.a){1+rc*exp[-2i*kz(z+hbeam)]}, but the terms in parentheses (rc and exp(...)) are
						 * replaced by their differences with -1 and 1, respectively (eliminating 1 inside)
						 */
						t2=ctemp*(rcp1*(t1+1)-t1);
						cvMultScal_RVec(t2,ex,b+j);
						if (which==INCPOL_X) {
							/* here we add (eIncRefl-ex*rc)*exp[ik(prIncRefl.r-2az*hbeam)] (non-zero only for
							 * p-polarization), expressing rc*exp(...) as t2-ctemp
							 */
							t3=2*(t2-ctemp);
							b[j]-=t3*ex[0];
							b[j+1]-=t3*ex[1];
						}
					}
				}
			}
			else for (i=0;i<local_nvoid_Ndip;i++) { // standard (non-surface) plane wave
				j=3*i;
				// can be replaced by complex multiplication (by precomputed phase factor), does not seem beneficial
				vSubtr(DipoleCoord+j,beam_center,r1);
				ctemp=imExp(WaveNum*DotProd(r1,prop)); // ctemp=exp(ik*r.a)
				cvMultScal_RVec(ctemp,ex,b+j); // b[i]=ctemp*ex
			}
			return;
		case B_DIPOLE: {
			double dip_p[3]; // dipole moment, = p0*prop
			vMultScal(p0,prop,dip_p);
			for (i=0;i<local_nvoid_Ndip;i++) { // here we explicitly use that dip_p is real
				j=3*i;
				vSubtr(DipoleCoord+j,beam_center,r1);
				(*InterTerm_real)(r1,gt);
				cSymMatrVecReal(gt,dip_p,b+j);
				if (surface) { // add reflected field
					r1[2]=DipoleCoord[j+2]+beam_center[2]+2*hsub;
					(*ReflTerm_real)(r1,gt);
					cReflMatrVecReal(gt,dip_p,v1);
					cvAdd(v1,b+j,b+j);
				}
			}
			/* calculate dipole inherent cross sections (for decay rate enhancements)
			 * General formula is C0=4pi*k*Im(p0*.G(r0,r0).p0) and it is used for reflection; but for direct
			 * interaction it is expected to be singular, so we use an explicit formula for point dipole:
			 * C0=(8pi/3)*k^4*|p0|^2. Thus we discard the choice of "-int ...", but still the used value should be
			 * correct for most formulations, e.g. poi, fcd, fcd_st, igt_so. Moreover, it is also logical, since the
			 * exciting dipole is really point one, in contrast to the dipoles composing the particle.
			 */
			double temp=p0*WaveNum*WaveNum; // in principle, t1=1/k, but we keep a general formula
			C0dipole=2*FOUR_PI_OVER_THREE*temp*temp;
			if (surface) {
				r1[0]=r1[1]=0;
				r1[2]=2*(beam_center[2]+hsub);
				(*ReflTerm_real)(r1,gt);
				double tmp;
				/* the following expression uses that dip_p is real and a specific (anti-)symmetry of the gt
				 * a general expression is commented out below
				 */
				tmp=dip_p[0]*dip_p[0]*cimag(gt[0])+2*dip_p[0]*dip_p[1]*cimag(gt[1])+dip_p[1]*dip_p[1]*cimag(gt[3])
				   +dip_p[2]*dip_p[2]*cimag(gt[5]);
//				v1[0]=dip_p[0]; v1[1]=dip_p[1]; v1[2]=dip_p[2];
//				cReflMatrVec(gt,v1,v2);
//				tmp=cDotProd_Im(v2,v1);
				C0dipole_refl=FOUR_PI*WaveNum*tmp;
			}
			return;
		}
		case B_LAGUERRE:
			/*
			. This 'GenerateB()' is called once (looped over every dipole point)
			. Returned at the end ONLY, where b is the E field returned (by reference, not directly)

			i 		=> each dipole index
			j=3*i 	=> jth vector start point
			ctemp   => value of E of beam at this point (store here temporarily) ---> doublecomplex type
						--> This is the magnitude of the E field essentially
			r1 		=> vector distance (double[3]) of dipole from origin (beam_center)
			DorProd(a,b) => dot product between two vectors (k=prop)
			DipoleCoord  => array of dipole coordinates;
			HENCE, DipoleCoord+j => coord of ith dipole -> [doublecomplex, doublecomplex, doublecomplex]
			cvMultScal_RVec()    => Multiplies a !doublecomplex SCALAR! by a !double VECTOR!, and stores the result in a !doublecomplex VECTOR!
			b => the output LIST of VECTORS for each dipole's E field
			ex, ey, prop 		 => The X,Y,Z vectors for the electric field (assuming prop is always in Z direction, true by default and in the cases important here)
			*/

			//printf("LOG PRINT; ---- Gen beam: LAGUERRE Overall ---- \n");
			double k, wavelength, w0, rho, z_R;
			double integral_start, integral_end;

			printf("local_nvoid_Ndip= %d \n",local_nvoid_Ndip);

			for (i=0;i<local_nvoid_Ndip;i++) { // standard (non-surface) plane wave
				/*
				l, then p specified in beam args
				
				- Bessel, Legendre, C all match
				- Python used alpha = 0.4*pow(10,-8)*(1.0+1j), 
							  beta  = 0.4*pow(10,-8)*(1.0-1j)
				*/

				//######################################
				//## MAYBE NEED TO SET VORTICITY HERE ##
				//######################################
				j=3*i;
				vSubtr(DipoleCoord+j,beam_center,r1);
				x=DotProd(r1,ex);		//Cartesian coords
				y=DotProd(r1,ey);		//
				z=DotProd(r1,prop);			//Cylindrical coords
				phi=atan2(y,x);				//
				rho = sqrt( (x*x) +(y*y) );	//

				//##########
				//### WANT TO PULL THE LAMBDA FROM THE DATA TO DO PROPERLY
				//##########
				k = WaveNum;	//In micrometres, is a double, NOT doublecomplex
				wavelength = (2.0*PI)/(k);
				w0 = 0.6*wavelength;			//<--- Was defined elsewhere for other beams (as a param, maybe worth doing same here)
				integral_start = 0.0;
				integral_end   = k;

				t1  = 1.0+1.0*I;//0.4*cpow(10,-8)*(1.0+1*I); 		//alpha
				t2  = 1.0-1.0*I;//0.4*cpow(10,-8)*(1.0-1*I); 		//beta
				z_R = k*(w0*w0)/2.0;	//z_R

				if(i==0) {
					printf("wavelength= %e \n", wavelength);
					printf("k = %e \n", k);
					printf("w0= %e \n", w0);
					printf("integral_start= %e \n", integral_start);
					printf("integral_end  = %e \n", integral_end);
					printf("t1 = %e, %e \n", creal(t1), cimag(t1));
					printf("t2 = %e, %e \n", creal(t2), cimag(t2));
					printf("z_R= %e \n", z_R);
				}

				//printf("bessel func jn(0,0)= %e \n",_jn(0));
				//printf("bessel func jn(0,1)= %e \n",_jn(1));

				t6 = integrate_simpsons_rule(integral_start, integral_end, 31, laguerre_getParam_E_component, 0, rho, phi, z, k, l, p, z_R, w0, t1, t2);	//E_x
				t7 = integrate_simpsons_rule(integral_start, integral_end, 31, laguerre_getParam_E_component, 1, rho, phi, z, k, l, p, z_R, w0, t1, t2);	//E_y
				t8 = integrate_simpsons_rule(integral_start, integral_end, 31, laguerre_getParam_E_component, 2, rho, phi, z, k, l, p, z_R, w0, t1, t2);	//E_z
				//--REAL PART--
				//t6 = creal(t6);//sqrt(creal(t6)*cimag(t6));	//Just take magnitude for E
				//t7 = creal(t7);//sqrt(creal(t7)*cimag(t7));	//
				//t8 = creal(t8);//sqrt(creal(t8)*cimag(t8));	//
				//--IMAGINARY PART--
				//t6 = cimag(t6);
				//t7 = cimag(t7);
				//t8 = cimag(t8);

				if(i==0) {
					printf("t6 = %e, %e \n", creal(t6), cimag(t6));
					printf("t7 = %e, %e \n", creal(t7), cimag(t7));
					printf("t8 = %e, %e \n", creal(t8), cimag(t8));
				}

				cvMultScal_RVec(t6,ex,v1);		//X component of E
				cvMultScal_RVec(t7,ey,v2);		//Y component of E
				cvMultScal_RVec(t8,prop,v3);	//Z component of E
				cvAdd2Self(v1,v2,v3);			//Bring together into v1
				ctemp = 1.0 +0.0*I;
				cvMultScal_cmplx(ctemp,v1,b+j);	//Setting E for dipole -> multiplied by 1.0 as a work around to set using cvMultScal_cmplx
			}
			return;
		case B_LMINUS:
		case B_DAVIS3:
		case B_BARTON5:
			for (i=0;i<local_nvoid_Ndip;i++) {
				j=3*i;
				// set relative coordinates (in beam's coordinate system)
				vSubtr(DipoleCoord+j,beam_center,r1);
				x=DotProd(r1,ex)*scale_x;
				y=DotProd(r1,ey)*scale_x;
				z=DotProd(r1,prop)*scale_z;
				ro2=x*x+y*y;
				Q=1/(2*z-I);
				psi0=-I*Q*cexp(I*Q*ro2);
				// ctemp=exp(ik*z0)*psi0, z0 - non-scaled coordinate (z/scale_z)
				ctemp=imExp(WaveNum*z/scale_z)*psi0;
				// the following logic (if-else-if...) is hard to replace by a simple switch
				if (beamtype==B_LMINUS) cvMultScal_RVec(ctemp,ex,b+j); // b[i]=ctemp*ex
				else {
					/* It is possible to rewrite the formulae below to avoid division by ro2, but we prefer
					 * dimensionless variables. The value for ro2=0 doesn't really matter (cancels afterwards).
					 * The current code should work OK even for very small ro2
					 */
					if (ro2==0) x2_s=0;
					else x2_s=x*x/ro2;
					Q2=Q*Q;
					ro4=ro2*ro2;
					// some combinations that are used more than once
					t4=s2*ro2*Q2; // t4=(s*ro*Q)^2
					t5=I*ro2*Q;   // t5=i*Q*ro^2
					t6=ro4*Q2;    // t6=ro^4*Q^2
					t7=x*s*Q;     // t7=x*s*Q
					if (beamtype==B_DAVIS3) {
						// t1=1+s^2(-4Q^2*x^2-iQ^3*ro^4)=1-t4(4x2_s+t5)
						t1 = 1 - t4*(4*x2_s+t5);
						// t2=0
						t2=0;
						// t3=-s(2Qx)+s^3(8Q^3*ro^2*x+2iQ^4*ro^4*x-4iQ^2x)=2t7[-1+iQ*s2*(-4t5+t6-2)]
						t3 = 2*t7*(-1 + I*Q*s2*(-4*t5+t6-2));
					}
					else if (beamtype==B_BARTON5) {
						if (ro2==0) xy_s=0; // see comment for x2_s above
						else xy_s=x*y/ro2;
						t8=8+2*t5; // t8=8+2i*Q*ro^2
						/* t1 = 1 + s^2(-ro^2*Q^2-i*ro^4*Q^3-2Q^2*x^2)
						 *    + s^4[2ro^4*Q^4+3iro^6*Q^5-0.5ro^8*Q^6+x^2(8ro^2*Q^4+2iro^4*Q^5)]
						 *    = 1 + t4*{-1-2x2_s-t5+t4*[2+3t5-0.5t6+x2_s*t8]}
						 */
						t1 = 1 + t4*(-1 - 2*x2_s - t5 + t4*(2+3*t5-0.5*t6+x2_s*t8));
						// t2=s^2(-2Q^2*xy)+s^4[xy(8ro^2*Q^4+2iro^4*Q^5)]=xy_s*t4(-2+t4*t8)
						t2=xy_s*t4*(-2+t4*t8);
						/* t3 = s(-2Qx) + s^3[(6ro^2*Q^3+2iro^4*Q^4)x] + s^5[(-20ro^4*Q^5-10iro^6*Q^6+ro^8*Q^7)x]
						 *    = t7{-2+t4[6+2t5+t4(-20-10t5+t6)]}
						 */
						t3 = t7*(-2 + t4*(6 + 2*t5 + t4*(-20-10*t5+t6)));
					}
					else LogError(ONE_POS,"Inconsistency in beam definition"); // to remove warnings
					// b[i]=ctemp(ex*t1+ey*t2+ez*t3)
					cvMultScal_RVec(t1,ex,v1);
					cvMultScal_RVec(t2,ey,v2);
					cvMultScal_RVec(t3,prop,v3);
					cvAdd2Self(v1,v2,v3);
					cvMultScal_cmplx(ctemp,v1,b+j);
				}
			}
			return;
#ifndef NO_FORTRAN
		case B_BES_CS:
		case B_BES_CSp:
		case B_BES_M:
		case B_BES_LE:
		case B_BES_LM:
		case B_BES_TEL:
		case B_BES_TML:
			/* we assume that matrix M determines x-polarization (in the beam reference frame), while y-one is obtained
			 * by rotation with additional phase shift
			 */
			vort=(which==INCPOL_Y) ? cpow(I,besN) : 1;
			for (i=0;i<local_nvoid_Ndip;i++) {
				j=3*i;
				vSubtr(DipoleCoord+j,beam_center,r1);
				x=DotProd(r1,ex);
				y=DotProd(r1,ey);
				z=DotProd(r1,prop);
				phi=atan2(y,x); // angular coordinate in a cylindrical coordinate system
				ctemp=imExp(besN*phi)*cexp(I*besKz*z)*vort/(WaveNum*WaveNum); // common factor
				arg=besKt*sqrt(x*x+y*y); // argument of Bessel functions
				if (arg<ROUND_ERR) {
					// TODO: the following seems incorrect for n=+-1 or +-2 (other fn will be non-zero)
					if (besN==0) fn[2]=1.;
					else fn[2]=0;
					fn[0]=fn[1]=fn[3]=fn[4]=0;
				}
				else {
					// TODO: the following looks very complicated, try to simplify (exp factors can be precalculated)
					n1=abs(besN)+2;
					if (besN<=-3) {
						bjndd_(&n1,&arg,jn1,td1,td2);
						if (n1%2==0) q=1;
						else q=-1;
						fn[0] =  q*jn1[-besN+2]*imExp(-2*phi);
						fn[1] = -q*jn1[-besN+1]*imExp(-phi);
						fn[2] =  q*jn1[-besN];
						fn[3] = -q*jn1[-besN-1]*imExp(phi);
						fn[4] =  q*jn1[-besN-2]*imExp(2*phi);
					}
					if (besN >= 2) {
						bjndd_(&n1,&arg,jn1,td1,td2);
						fn[0] = jn1[besN-2]*imExp(-2*phi);
						fn[1] = jn1[besN-1]*imExp(-phi);
						fn[2] = jn1[besN];
						fn[3] = jn1[besN+1]*imExp(phi);
						fn[4] = jn1[besN+2]*imExp(2*phi);
					}
					if (besN == -2) {
						bjndd_(&n1,&arg,jn1,td1,td2);
						fn[0] =  jn1[4]*imExp(-2*phi);
						fn[1] = -jn1[3]*imExp(-phi);
						fn[2] =  jn1[2];
						fn[3] = -jn1[1]*imExp(phi);
						fn[4] =  jn1[0]*imExp(2*phi);
					}
					if (besN == -1) {
						bjndd_(&n1,&arg,jn1,td1,td2);
						fn[0] = -jn1[3]*imExp(-2*phi);
						fn[1] =  jn1[2]*imExp(-phi);
						fn[2] = -jn1[1];
						fn[3] =  jn1[0]*imExp(phi);
						fn[4] =  jn1[1]*imExp(2*phi);
					}
					if (besN == 0) {
						bjndd_(&n1,&arg,jn1,td1,td2);
						fn[0] =  jn1[2]*imExp(-2*phi);
						fn[1] = -jn1[1]*imExp(-phi);
						fn[2] =  jn1[0];
						fn[3] =  jn1[1]*imExp(phi);
						fn[4] =  jn1[2]*imExp(2*phi);
					}
					if (besN == 1) {
						bjndd_(&n1,&arg,jn1,td1,td2);
						fn[0] = -jn1[1]*imExp(-2*phi);
						fn[1] =  jn1[0]*imExp(-phi);
						fn[2] =  jn1[1];
						fn[3] =  jn1[2]*imExp(phi);
						fn[4] =  jn1[3]*imExp(2*phi);
					}
				}
				/* TODO: all factors before f[n] should be calculated in InitBeam beforehand
				 * this will improve speed and allow robust calculation in the limit, e.g., of kt->0
				 */
				t1 = (((WaveNum*WaveNum+besKz*besKz)/2.*besM[0] +  WaveNum*besKz*besM[3])*fn[2] +
					  besKt*besKt/4.*((besM[0]+I*besM[1])*fn[0] + (besM[0]-I*besM[1])*fn[4])); // Ex
				t2 = (((WaveNum*WaveNum+besKz*besKz)/2.*besM[1] -  WaveNum*besKz*besM[2])*fn[2] +
					  I*besKt*besKt/4.*((besM[0]+I*besM[1])*fn[0] - (besM[0]-I*besM[1])*fn[4])); // Ey
				t3 = ((I*besKz*(besM[0]+I*besM[1]) + WaveNum*(besM[2]+I*besM[3]))*fn[1] -
					  (I*besKz*(besM[0]-I*besM[1]) - WaveNum*(besM[2]-I*besM[3]))*fn[3])*besKt/2.; // Ez

				cvMultScal_RVec(t1,ex,v1);
				cvMultScal_RVec(t2,ey,v2);
				cvMultScal_RVec(t3,prop,v3);
				cvAdd2Self(v1,v2,v3);
				cvMultScal_cmplx(ctemp,v1,b+j);
			}
			return;
#endif // !NO_FORTRAN
		case B_READ:
			if (which==INCPOL_Y) fname=beam_fnameY;
			else fname=beam_fnameX; // which==INCPOL_X
			ReadField(fname,b);
			return;
	}
	LogError(ONE_POS,"Unknown type of incident beam (%d)",(int)beamtype);
	/* TO ADD NEW BEAM
	 * add a case above. Identifier ('B_...') should be defined inside 'enum beam' in const.h. This case should set
	 * complex vector 'b', describing the incident field in the particle reference frame. It is set inside the cycle for
	 * each dipole of the particle and is calculated using
	 * 1) 'DipoleCoord' - array of dipole coordinates;
	 * 2) 'prop' - propagation direction of the incident field;
	 * 3) 'ex' - direction of incident polarization;
	 * 4) 'ey' - complementary unity vector of polarization (orthogonal to both 'prop' and 'ex');
	 * 5) 'beam_center' - beam center in the particle reference frame (automatically calculated from 'beam_center_0'
	 *                    defined by '-beam_center' command line option).
	 * If the new beam type is compatible with '-surf', include here the corresponding code. For that you will need
	 * the variables, related to surface - see vars.c after "// related to a nearby surface".
	 * If you need temporary local variables (which are used only in this part of the code), either use 't1'-'t8' or
	 * define your own (with more informative names) in the beginning of this function.
	 */
}