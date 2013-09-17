StellarSpiral
=============

The project contains excution files and library files.


## Excution Files 

* **Corotation.exe**
Plot Omega, Omega+-Kappa graph. 
Show corotation,pattern speed and 4-1 resonance.
		
* **Density.exe**
Draw stellar density distribution.

		usage:   Density.exe [option]            
		options:                         
		-p, --project [degree]  To project the density.
		-c, --circle            To draw circles.       
		-r, --zauto             Automatic find scale of z
		-z, --zmax [scale of z] Specified the scale of z
		-h, --help              Show this help page     

* **Ecg.exe**
Show amplitude of some radius.

		usage:  Ecg.exe -i filename.h5 [-z zmax]
         		Ecg.exe -v
         		
* **FindOne.exe**
Find pattern speed. 

		usage: FindOne.exe wr wi
		
* **Follow.exe**
Measure amplitute along two spirals. 

		usage:   Follow.exe [options] -i Filename.h5
		options:
         -z, --zmax [scale of z] Specified the scale of z
         -h, --help              Show this help page
         -v,                     Version information.
         -x,                     Number of XWindows.
         -c,                     Plot circles 3,4,5.
         -s,                     Plot stellar density.
         --output                Save found points.

* **GasDensity.exe**
Draw gas density distribution.

		usage:   GasDensity.exe [option]
		options:
         -c, --circle            To draw circles.
         -i [file name.h5]       Read file.
         -m [start] [end]        Read files.
         -p, --project [degree]  To project the density.
         -r, --zauto             Automatic find scale of z
         -z, --zmax [scale of z] Specified the scale of z
         -h, --help              Show this help page
         -s, --stellar           Draw stellar contour.
         -q,                     Draw instability.
         -v,                     Version information.
         --save                  Save as png file.
         --cut                   Make a cut at CO.


* **Kendall.exe**
Calculate prefactor of density contrast.

		usage:   Kendall.exe [option]
		options:
         -f, --force             Increase force multiple.

* **SearchAll.exe**
Find initial guess for FindOne.exe
Using **SearchAllMpi.sh** to submit to dl clusters.

* **Shift.exe**
Trace the amplitude of both stellar and gaseous spirals.

		usage:   Follow.exe [options] -i Filename.h5
		options:
         -z, --zmax [scale of z] Specified the scale of z
         -h, --help              Show this help page
         -v,                     Version information.
         -x,                     Number of XWindows.
         -c,                     Plot circles at 3,4 and 5 kpc.
         --output                Output Phase Shift.

## note
Use Kendal.exe to find the amplitude factor.
The amp parameter in antares/para.list and the factor should be the same.
