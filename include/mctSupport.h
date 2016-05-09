#ifndef _INC_MCT_SUPPORT_H
#define _INC_MCT_SUPPORT_H

#include "mctGlobals.h"

namespace mct
{

inline double convertDegreesToRadians (double x)
{ 
	return (x * (PI/180.)); 
}

inline double convertRadiansToDegrees (double x)
{ 
	return (x*(180./PI)); 
}

}
#endif