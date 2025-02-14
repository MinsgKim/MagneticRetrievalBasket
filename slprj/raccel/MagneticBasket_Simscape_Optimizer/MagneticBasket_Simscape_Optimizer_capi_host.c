#include "MagneticBasket_Simscape_Optimizer_capi_host.h"
static MagneticBasket_Simscape_Optimizer_host_DataMapInfo_T root;
static int initialized = 0;
__declspec( dllexport ) rtwCAPI_ModelMappingInfo *getRootMappingInfo()
{
    if (initialized == 0) {
        initialized = 1;
        MagneticBasket_Simscape_Optimizer_host_InitializeDataMapInfo(&(root), "MagneticBasket_Simscape_Optimizer");
    }
    return &root.mmi;
}

rtwCAPI_ModelMappingInfo *mexFunction(){return(getRootMappingInfo());}
