/* SpikeFilter.c - gets spikes and LFP from Intan Data File(s)
 *
 * The calling syntax is:
 *
 *      Sf = SpikeFilter(S, min, max, fs)
 *
 * Compile with:
 *
 *      mex SpikeFilter.c
 *
 * INPUT
 *  S - Nsamples-by-Nchannels array of LFP data to filter for spikes
 *  min - min freq cutoff (or min wavelet scale) to use for filtering
 *  max - max freq cutoff (or upper wavelet scale) to use for filtering
 *  fs - sampling frequency (use -1 for wavelet filter)
 *
 * OUTPUT
 *  Sf - spike-filtered S
 *
 * This is a MEX file for MATLAB.
*/

#include "mex.h"
#include <math.h>

const double M_PI = 3.1415962;
void wtfilt(char, unsigned int, double **, double **); //filter generator


/* Arbitrary Partial Discrete Wavelet Transform
 *
 * INPUT
 *	a - the vector to be transformed
 *	n - length of a
 *	isign - if >0 will do forward wavelet transform, else inverse WT
 *	nc - number of filter values
 *	Ds - smoothing filter values
 *	Dd - detail filter values
 *
 * OUTPUT
 *	replaces a[] with its wavelet transform
 */
void pwt(double a[], unsigned int n, int isign, unsigned int nc, double * Ds, double * Dd)
{
    // Declare variables
    double ai, ai1, * wksp_zo, * wksp;
	unsigned int i,ii,jf,k,ni,nh;
    int off;
    
    // Base case
    if (n < 4) return;
    
    // Initialize constants
    wksp_zo = malloc(n*sizeof(double)); //to store output
    wksp = wksp_zo - 1; //unit-offset pointer to output array
    off = -( (int) (nc >> 1) ); //filter offset
    nh = n>>1; //half of n
    for (i=1; i<=n; i++) wksp[i]=0.0; //initialize outupt to 0
    
    // Forward or inverse transform
    if (isign >= 0) { //forward transform
        for (ii=1,i=1; i<=n; i+=2, ii++) {
            ni = i+off+n;
            for (k=1; k<=nc; k++) { //calculate coeff @ this location
                jf = (ni+k)%n+1;
                wksp[ii] += Ds[k]*a[jf];
                wksp[ii+nh] += Dd[k]*a[jf];
            }
        }
    } else { //inverse transform
        for (ii=1,i=1; i<=n; i+=2,ii++) {
            ai = a[ii];
            ai1 = a[ii+nh];
            ni = i+off+n;
            for (k=1; k<=nc; k++) { //calculate signal value @ this location
                jf = (ni+k)%n+1;
                wksp[jf] += Ds[k]*ai + Dd[k]*ai1;
            }
        }
    }
    
    // Write to output array
    for (i=1; i<=n; i++) a[i] = wksp[i];
    free(wksp_zo);
}


/* Does either forward or inverse discrete wavelet transform using 
 *  pyramidal scheme and the specified wavelet.
 *
 * INPUT
 *  a - vector to be transformed
 *  n - length of a[]
 *  isign - if isign>0 do forward transform, else to inverse transform
 *  nc - number of wavelet filter coefficients
 *  wtype - a single character indicating the type of filter to use
 *      Available types:
 *          d - Daubechies wavelet
 *          c - Coiflet wavelet
 * 
 * OUTPUT
 *  Replaces a[] with its wavelet transform
 */
void dwt(double a[], unsigned int n, int isign, unsigned int nc, char wtype)
{
    // Declare variables
	unsigned int nn;
	double * Ds, *Dd; //pointer to array of filter coefficients

    // Base case
    if (n<nc) return;

    // Get the filter coefficients for desired filter
    wtfilt(wtype, nc, &Ds, &Dd);

    // Do the transform using the pyramidal scheme
    if (isign >= 0) { //forward wavelet transform
        for (nn=n; nn>=4; nn>>=1) pwt(a-1, nn, isign, nc, Ds-1, Dd-1);
    } else { //inverse wavelet transform
        for (nn=4; nn<=n; nn<<=1) pwt(a-1, nn, isign, nc, Ds-1, Dd-1);
    }

    // Free the filter vectors
    free(Ds); 
    free(Dd);

}


/* RcBandpassFilt
    
Bandpass a signal using an RC filter

INUPTS
    S - the original unfiltered signal
    Sf - the output filtered signal
    N - number of samples in S
    lc - lower cutoff frequency (Hz)
    hc - upper cutoff frequency (Hz)
    fs - sampling frequency (Hz)
 *
OUTPUTS
    Fills Sf with a highpassed version of S
*/
void RcBandpassFilt(double * S, double * Sf, 
        unsigned int N, double lc, double hc, double fs)
{
    // Declare Variables
    int i;
    double RC, a, a1, dt;
    
    // Highpass
    dt = 1.0/fs;
    RC = 1.0/(2.0*M_PI*lc);
    a = RC / (RC + dt);
    Sf[0] = S[0];
    for (i=1; i<N; i++) {
        Sf[i] = a*(Sf[i-1] + S[i] - S[i-1]);
    }
    
    // Lowpass
    RC = 1.0/(2.0*M_PI*hc);
    a = dt / (RC + dt); //alpha
    a1 = 1.0-a;
    for (i=1; i<N; i++) {
        Sf[i] = a*Sf[i] + a1*Sf[i-1];
    }
        
}


/* MATLAB MEX FUNCTION: Sf = SpikeFilter(S, min, max, fs) */
void mexFunction( int nlhs, mxArray *plhs[], 
                  int nrhs, const mxArray * prhs[])
{
    // Declare variables and array pointers
    double * S;
    double * S2;
    const mwSize * dims;
    unsigned int Ns, Nc;
    double lowcut, highcut, fs;
    double * Sf;
    int iC, i;
    unsigned int numC, np2;
    unsigned int maxS, minS;
    
    
    // Get input from matlab
    S = mxGetPr(prhs[0]); //raw LFP
    dims = mxGetDimensions(prhs[0]);
    Ns = (unsigned int) dims[0];  //number of (data) points per channel
    Nc = (unsigned int) dims[1];  //number of channels
    fs = (double) mxGetScalar(prhs[3]); //sampling frequency
    
    // Allocate output array
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    Sf = mxGetPr(plhs[0]);
    
    // RC filter
    if (fs>0) {
        lowcut = (double) mxGetScalar(prhs[1]); //for RC bandpass
        highcut = (double) mxGetScalar(prhs[2]);
        for (iC=0; iC<Nc; iC++) {
            RcBandpassFilt(S+Ns*iC, Sf+Ns*iC, Ns, lowcut, highcut, fs);
        }
        
    } else { // Wavelet filter
        minS = (unsigned int) mxGetScalar(prhs[1]); //scales for 
        maxS = (unsigned int) mxGetScalar(prhs[2]); //wavelet filter
        np2 = pow(2, ceil(log(Ns)/log(2.0))); //next power of 2 from Ns
        //NOTE: the pow and log functions take doubles!!!
        S2 = malloc(np2*sizeof(double)); //next-power-of-2 length array
        numC = 18; //number of wavelet coeffs to use (can be 18, 12, or 6 for coiflet)
        for (iC=0; iC<Nc; iC++) {
            for (i=0; i<Ns; i++) S2[i] = S[Ns*iC+i]; //copy S
            for (i=Ns; i<np2; i++) S2[i] = 0; //pad with 0s
            dwt(S2, np2, 1, numC, 'c'); //forward transform
            for (i=0; i<(np2>>maxS); i++) S2[i] = 0; //remove unwanted scales
            for (i=(np2>>(minS-1)); i<np2; i++) S2[i] = 0;
            dwt(S2, np2, -1, numC, 'c'); //reconstruct signal
            for (i=0; i<Ns; i++) Sf[Ns*iC+i] = S2[i]; //copy to output
        }
    }
        
}


/* Generates wavelet filter coefficients 
 * 
 * INPUT
 *  wtype - character indicating the type of filter desired (see FILTER TYPES)
 *  nc - number of filter coefficients
 *  Fs - address of a double pointer, will set to point to smoothing filter coeffs
 *  Fd - address of a double pointer, will set to point to detail filter coeffs
 *
 * FILTER TYPES
 *	d - Daubechies wavelet
 *	c - Coiflet wavelet
 */
void wtfilt(char wtype, unsigned int nc, double ** Fs, double ** Fd)
{
    // Declare variables
    double * Ds, * Dd, sig;
    unsigned int i;
    
    // Allocate filter vectors
    Ds = malloc(nc*sizeof(double));
    Dd = malloc(nc*sizeof(double));
    
    // Define the coefficients
    switch ( wtype ) { //what filter to use?
        case 'd' : //Daubechies wavelets
            switch ( nc ) { 
                case 2 :  // Daubechies 2-coeff (aka Haar wavelet)
                    Ds[0] = 1;
                    Ds[1] = 1;
                    break;
                case 4 :  // Daubechies 4-coeff (are these backwards?)
                    Ds[0] = 0.4829629131445341;
                    Ds[1] = 0.8365163037378079;
                    Ds[2] = 0.2241438680420134;
                    Ds[3] = -0.1294095225512604;
                    break;
                case 6 : // Daubechies 6-coeff
                    Ds[0] = 0.33267055295008261;
                    Ds[1] = 0.80689150931109257;
                    Ds[2] = 0.45987750211849157;
                    Ds[3] = -0.13501102001025458;
                    Ds[4] = -0.08544127388202666;
                    Ds[5] = 0.03522629188570953;
                    break;
                case 8 : // Daubechies 8-coeff
                    Ds[0] = 0.2303778133088965008632911830440708500016152482483092977910968;
                    Ds[1] = 0.7148465705529156470899219552739926037076084010993081758450110;
                    Ds[2] = 0.6308807679298589078817163383006152202032229226771951174057473;
                    Ds[3] = -0.02798376941685985421141374718007538541198732022449175284003358;
                    Ds[4] = -0.1870348117190930840795706727890814195845441743745800912057770;
                    Ds[5] = 0.03084138183556076362721936253495905017031482172003403341821219;
                    Ds[6] = 0.03288301166688519973540751354924438866454194113754971259727278;
                    Ds[7] = -0.01059740178506903210488320852402722918109996490637641983484974;
                    break;
                default:
                    mexPrintf("ERROR: unimplemented number of filter coefficients for Daubechies wavelet\n");
                    return;        
            }
            break;
        case 'c' : //Coiflet wavelets
            switch ( nc ) { 
                case 6 :  // Coiflet with 1 vanishing moment
                    Ds[0] = -0.015655728;
                    Ds[1] = -0.07273262;
                    Ds[2] = 0.384864847;
                    Ds[3] = 0.85257202;
                    Ds[4] = 0.337897662;
                    Ds[5] = -0.07273262;
                    break;
                case 12 :  // Coiflet with 2 vanishing moments
                    Ds[0] = -0.000720549;
                    Ds[1] = -0.001823209;
                    Ds[2] = 0.005611435;
                    Ds[3] = 0.023680172;
                    Ds[4] = -0.059434419;
                    Ds[5] = -0.076488599;
                    Ds[6] = 0.417005184;
                    Ds[7] = 0.812723635;
                    Ds[8] = 0.386110067;
                    Ds[9] = -0.067372555;
                    Ds[10] = -0.041464937;
                    Ds[10] = 0.016387336;
                    break;
                case 18 :  // Coiflet with 3 vanishing moments
                    Ds[0] = -0.0000346;
                    Ds[1] = -0.000071;
                    Ds[2] = 0.000466217;
                    Ds[3] = 0.001117519;
                    Ds[4] = -0.002574518;
                    Ds[5] = -0.009007976;
                    Ds[6] = 0.015880545;
                    Ds[7] = 0.034555028;
                    Ds[8] = -0.082301927;
                    Ds[9] = -0.071799822;
                    Ds[10] = 0.428483476;
                    Ds[11] = 0.793777223;
                    Ds[12] = 0.405176902;
                    Ds[13] = -0.06112339;
                    Ds[14] = -0.065771911;
                    Ds[15] = 0.023452696;
                    Ds[16] = 0.007782596;
                    Ds[17] = -0.003793513;
                    break;
                default:
                    mexPrintf("ERROR: unimplemented number of filter coefficients for Coiflet wavelet\n");
                    return;
            }
            break;
        default: {
            mexPrintf("ERROR: unimplemented wavelet type: %c\n", wtype);
            return;
        }
    }
            
    // Define detail coefficients
    sig = -1.0;
    for (i=0; i<nc; i++) {
        Dd[nc-i-1] = sig * Ds[i];
        sig = -sig;
    }
    
    // Set the return pointers
    *(Fs) = Ds;
    *(Fd) = Dd;
}
