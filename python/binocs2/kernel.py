# BINOCS OpenCL kernel subroutines

def sedkernel(nopt, nnir, nmir, type="default"):
	"""SEDKERNEL
	DESCRIPTION: Returns chosen OpenCL SED kernel
	INPUT:       nopt -- number of good optical filters
	             nnir -- number of good near IR filters
	             nmir -- number of good mid IR filters
	             type -- (optional) string determining type of kernel to be returned
	                   "chi", "default"
	OUTPUT:      kernelstr -- string containing C program to be compiled
	"""
	
	# Chi^2 Kernel
	if type == "chi":
		print("Using Chi^2 Kernel")
		kernelstr = """
			__kernel void binsub( __global float* iso, __global float* data, __global float* err, __global float* chi, __global int* fit, const float chithresh, const int nmodels ) {
				int s = get_global_id(0);
				fit[s] = -1.0;
				chi[s] = -1.0;
				float bestchi = 1000.0;
				int bestfit = -1;
		
				// Loop through models
				for (int m = 0; m < nmodels; m++){
					// Initialize variables for this run
					float tmpchi = 0.0, thischi = 0.0, totfilt = 0.0;
					int gubv = 0, gsds = 0, gvis = 0, gnir = 0, gmir = 0;
	
					// Loop through filters and compare star to the model
					for (int f = 0; f < 17; f++){
						thischi = ((data[17*s+f] - iso[23*m+f+6]) * (data[17*s+f] - iso[23*m+f+6])) / (err[17*s+f] * err[17*s+f]);
						if (thischi < chithresh){
							if (f < 5){ gubv++; }
							else if (f < 10){ gsds++; }
							else if (f < 13){ gnir++; }
							else{ gmir++; }
							totfilt++;
							tmpchi += thischi;
						} 
						// If star is more than 100x the uncertainty away on this filter it *will not* fit the star. Abort.
						else if (thischi > 100 && data[17*s+f] < 80) { break; }
					}
					// See which visual filter set has more matches
					if (gubv > gsds){ gvis = gubv; }
					else {gvis = gsds; }
					// See if this comparison has enough filters to be used
					if (gvis >= %d && gnir >= %d && gmir >= %d){
						// See if this model is better than the previous best
						if (tmpchi / totfilt < bestchi){
							bestchi = tmpchi / totfilt;
							bestfit = m;
						}
					}
				}
				// Save best-fit model
				chi[s] = bestchi;
				fit[s] = bestfit;
			}""" % (nopt, nnir, nmir)
	
	# Default isochrone comparison kernel
	else: 
		kernelstr = """
			__kernel void binsub( __global float* iso, __global float* data, __global float* err, __global float* chi, __global int* fit, const float chithresh, const int nmodels ) {
				int s = get_global_id(0);
				fit[s] = -1.0;
				chi[s] = -1.0;
				float bestchi = -1.0;
				int bestfit = -1;
		
				// Loop through models
				for (int m = 0; m < nmodels; m++){
					// Initialize variables for this run
					float tmpchi = 0.0, thischi = 0.0, totfilt = 0.0;
					int gubv = 0, gsds = 0, gvis = 0, gnir = 0, gmir = 0;
			
					// Loop through filters and compare star to the model
					for (int f = 0; f < 17; f++){
						thischi = 1.0 / (fabs(data[17*s+f] - iso[23*m+f+6]) + 0.01);
						if (thischi > chithresh){
							if (f < 5){ gubv++; }
							else if (f < 10){ gsds++; }
							else if (f < 13){ gnir++; }
							else{ gmir++; }
							totfilt++;
							tmpchi += thischi;
						} 
						// If star is more than 2 magnitudes away on this filter it *will not* fit the star. Abort.
						else if (thischi < 0.5 && data[17*s+f] < 80) { break; }
					}
					// See which visual filter set has more matches
					if (gubv > gsds){ gvis = gubv; }
					else {gvis = gsds; }
					// See if this comparison has enough filters to be used
					if (gvis >= %d && gnir >= %d && gmir >= %d){
						// See if this model is better than the previous best
						if (tmpchi / totfilt > bestchi){
							bestchi = tmpchi / totfilt;
							bestfit = m;
						}
					}
				}
				// Save best-fit model
				chi[s] = bestchi;
				fit[s] = bestfit;
			}""" % (nopt, nnir, nmir)
			
	return kernelstr