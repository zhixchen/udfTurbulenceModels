# udfTurbulenceModels
## Prerequisites
* OpenFOAM version 7 

## Getting started:

  * load the OpenFOAM environment.
   
  * Download the code and compile using 'wmake'.
  
  * Add the following line into your 'controlDict':
  
    `libs ("libudfTurbulenceModels.so")`
    
  * In 'constant/turbulenceProperties', e.g. for Sigma model, set:
  
    `LESModel    Sigma`
    
    `SigmaCoeffs{ filter  simple; }`
    
