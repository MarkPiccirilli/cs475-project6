/********************************
* Title: Cuda Monte Carlo
* Name: Mark Piccirilli
* Course: CS 475
* Assignment: project 6
* Last Modified:
* Description:
* ************************************/

// System includes
#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>

// CUDA runtime
#include <cuda_runtime.h>

// Helper functions and utilities to work with CUDA
#include "helper_functions.h"
#include "helper_cuda.h"


#ifndef BLOCKSIZE
#define BLOCKSIZE		32		// number of threads per block
#endif

#ifndef NUMTRIALS
#define NUMTRIALS		64000		// to make the timing more accurate
#endif

#ifndef NUMTRIES
#define NUMTRIES		100		// to make the timing more accurate
#endif

#ifndef TOLERANCE
#define TOLERANCE		0.00001f	// tolerance to relative error
#endif

//ranges for random variables
const float XCMIN = 0.0;
const float XCMAX = 2.0;
const float YCMIN = 0.0;
const float YCMAX = 2.0;
const float RMIN = 0.5;
const float RMAX = 2.0;

// Monte Carlo Simulation

__global__  void MonteCarlo( float *xcs, float *ycs, float *rs, float* numHits )
{
	__shared__ float hits[BLOCKSIZE];

	unsigned int numItems = blockDim.x;
	unsigned int tnum = threadIdx.x;
	unsigned int wgNum = blockIdx.x;
	unsigned int gid = blockIdx.x*blockDim.x + threadIdx.x;

	//start out assuming its a hit.  If it is not a hit, adjust hits[tnum] to 0
	hits[tnum] = 1;

	float xc = xcs[gid];
	float yc = ycs[gid];
	float a = 2.;
	float b = -2.*(xc + yc);
	float c = xc*xc + yc*yc - r*r;
	float d = b*b - 4.*a*c;

	if(d < 0) {
		//laser misses circle
		hits[tnum] = 0;
	}

	d = sqrt(d);
	float t1 = (-b + d)/(2.*a);
	float t2 = (-b - d)/(2.*a);
	float tmin;
	//float tmin = t1 < t2 ? t1 : t2;
	if(t1 < t2) {
		tmin = t1;
	} 
	else {
		tmin = t2;
	}

	if(tmin < 0) {
		//circle engulfs laser
		hits[tnum] = 0;
	}

	float xcir = tmin;
	float ycir = tmin;

	//unitize normal vector
	float normalx = xcir - xc;
	float normaly = ycir - yc;
	float normal = sqrt(normalx*normalx + normaly*normaly);
	//unit vectors
	normalx /= normal;
	normaly /= normal;

	//unitized incoming vector
	float inx = xcir - 0.;
	float iny = ycir - 0.;
	float in = sqrt(inx*inx + iny*iny);

	//bounced vector
	float dot = inx*normalx + iny*normaly;
	float outx = inx - 2.*normalx*dot;
	float outy = iny - 2.*normaly*dot;

	float t = (0. - ycir)/outy;

	if(t < 0.) {
		//beam went up
		hits[tnum] = 0;
	}

	//add together all hits in block using reduction
	for (int offset = 1; offset < numItems; offset *= 2)
	{
		int mask = 2 * offset - 1;
		__syncthreads();
		if ((tnum & mask) == 0)
		{
			hits[tnum] += hits[tnum + offset];
		}
	}

	__syncthreads();
	if (tnum == 0)
		numHits[wgNum] = hits[0];
}


// main program:

int
main( int argc, char* argv[ ] )
{
	int dev = findCudaDevice(argc, (const char **)argv);

	// allocate host memory:

	float * hxcs = new float [ NUMTRIALS ];
	float * hycs = new float [ NUMTRIALS ];
	float * hrc = new float [ NUMTRIALS ];
	int * hnumHits = new int [ NUMTRIALS/BLOCKSIZE ];

	for( int i = 0; i < NUMTRIALS; i++ )
	{
		hxcs[n] = Ranf(XCMIN, XCMAX);
		hycs[n] = Ranf(YCMIN, YCMAX);
		rs[n] = Ranf(RMIN, RMAX);
	}

	// allocate device memory:

	float *dxcs, *dycs, *drs, *dnumHits;

	dim3 dimsA( NUMTRIALS, 1, 1 );
	dim3 dimsB( NUMTRIALS, 1, 1 );
	dim3 dimsC( NUMTRIALS, 1, 1 );
	dim3 dimsD( SIZE/BLOCKSIZE, 1, 1 );

	//__shared__ float prods[SIZE/BLOCKSIZE];


	cudaError_t status;
	status = cudaMalloc( reinterpret_cast<void **>(&dxcs), NUMTRIALS*sizeof(float) );
		checkCudaErrors( status );
	status = cudaMalloc( reinterpret_cast<void **>(&dycs), NUMTRIALS*sizeof(float) );
		checkCudaErrors( status );
	status = cudaMalloc( reinterpret_cast<void **>(&drsx), NUMTRIALS*sizeof(float) );
		checkCudaErrors( status );
	status = cudaMalloc( reinterpret_cast<void **>(&dnumHits), (SIZE/BLOCKSIZE)*sizeof(int) );
		checkCudaErrors( status );


	// copy host memory to the device:

	status = cudaMemcpy( dxcs, hxcs, NUMTRAILS*sizeof(float), cudaMemcpyHostToDevice );
		checkCudaErrors( status );
	status = cudaMemcpy( dycs, hycs, NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );
		checkCudaErrors( status );
	status = cudaMemcpy( drs, hrs, NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );
		checkCudaErrors( status );

	// setup the execution parameters:

	dim3 threads(BLOCKSIZE, 1, 1 );
	dim3 grid( NUMTRIALS / threads.x, 1, 1 );

	// Create and start timer

	cudaDeviceSynchronize( );

	// allocate CUDA events that we'll use for timing:

	cudaEvent_t start, stop;
	status = cudaEventCreate( &start );
		checkCudaErrors( status );
	status = cudaEventCreate( &stop );
		checkCudaErrors( status );

	// record the start event:

	status = cudaEventRecord( start, NULL );
		checkCudaErrors( status );

	// execute the kernel:

	for( int t = 0; t < NUMTRIALS; t++)
	{
	        MonteCarlo<<< grid, threads >>>( dxcs, dycs, drcs, dnumHits );
	}

	// record the stop event:

	status = cudaEventRecord( stop, NULL );
		checkCudaErrors( status );

	// wait for the stop event to complete:

	status = cudaEventSynchronize( stop );
		checkCudaErrors( status );

	float msecTotal = 0.0f;
	status = cudaEventElapsedTime( &msecTotal, start, stop );
		checkCudaErrors( status );

	// compute and print the performance

	double secondsTotal = 0.001 * (double)msecTotal;
	double trailsPerSecond = (float)NUMTRIALS * (float)NUMTRIES / secondsTotal;
	double megaTrailsPerSecond = trailsPerSecond / 1000000.;
	fprintf( stderr, "Number of trials = %10d, MegaTrials/Second = %10.2lf\n", NUMTRIALS, megaTrialsPerSecond );

	// copy result from the device to the host:

	status = cudaMemcpy( hnumHits, dnumHits, (SIZE/BLOCKSIZE)*sizeof(int), cudaMemcpyDeviceToHost );
		checkCudaErrors( status );

	// check the sum :

	int numHits = 0;
	for(int i = 0; i < NUMTRIALS/BLOCKSIZE; i++ )
	{
		//fprintf(stderr, "hnumHits[%6d] = %d\n", i, hnumHits[i]);
		numHits += hnumHits[i];
	}
	fprintf( stderr, "\nnumHits = %ld\n", numHits );
	
	//calculate frequency
	float frequency = numHits/NUMTRIALS;

	printf("frequency = %lf\n", frequency);

	// clean up memory:
	delete [ ] hxcs;
	delete [ ] hycs;
	delete [ ] hrs;
	delete [ ] hnumHits;

	status = cudaFree( dxcs );
		checkCudaErrors( status );
	status = cudaFree( dycs );
		checkCudaErrors( status );
	status = cudaFree( drs );
		checkCudaErrors( status );
	status = cudaFree( dnumHits );
		checkCudaErrors( status );


	return 0;
}

