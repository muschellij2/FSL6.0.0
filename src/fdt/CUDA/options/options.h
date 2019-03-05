#define SIZE_SUB_PART 12800 	// Number of voxels to compute in each subpart of the volume. 12800 = 16 SM * 800. It uses about 500MB of global memory
				// Could increase it, but it takes almost the same time. 
				// Enough to exploit parallelism and low memory requirements for old GPUs

#define THREADS_VOXEL_MCMC 32   // Multiple of 32: Threads collaborating to compute a voxel
#define VOXELS_BLOCK_MCMC 8	// Number of voxels per Block

#define THREADS_BLOCK_FIT 64
#define THREADS_BLOCK_RAND 256


