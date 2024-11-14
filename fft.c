#include <math.h>
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

'''
nn is the length of the input (i.e., nn=N) which should be a power of 2.
Warning: the program does not check this condition).

The elements in the array data[] (i.e., f(x) for the forward FFT or F(u) for
the inverse FFT) must be stored from data[1] to data[2*nn]; data[0] is not
used.

The real part of the input is stored in the odd locations of the array
(data[1], data[3], data[5]. etc) and the imaginary part in the even locations
(data[2], data[4], data[6], etc.). Please note that the imaginary part is
zero when you compute the forward FFT of f(x). This not the case when you
compute the inverse FFT (i.e., the imaginary part of F(u) is typically not
zero). In general, both f(x) and F(u) can be thought as being complex
functions.

isign: -1 Forward FFT, isign: 1 Inverse FFT (based on our definition)
Warning: the FFT routine provided does not multiply by the normalization
factor 1/N that appears in the forward DFT equation; you should do this
yourself.
'''



void fft(float data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software 0#Y". */


int main() {
	std::cout << "test"


	return 0
}