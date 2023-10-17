// The way to end the program gracefully is different for each main program.
// In Standalone, the Graceful Ends depend on the global variable pModel to output the results to
//   the output file defined in the rvi file.
// In BMI, the Graceful Ends does not depend on the global variable pModel and output the results to
//   the standard output.

#ifdef STANDALONE
    #include "GracefulEndStandalone.h"
#elif BMI_LIBRARY
    #include "GracefulEndBMI.h"
#endif
