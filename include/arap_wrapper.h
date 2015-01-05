#ifndef ARAP_WRAPPER_H
#define ARAP_WRAPPER_H

#ifdef WIN32
//----------------------------------------------------------------------
//  To compile the code into a windows DLL, you must define the 
//  symbol ARAP_DLL_EXPORTS. 
// 
//  To compile the code statically into a windows executable 
//    (i.e. not using a separate DLL) define ARAP_DLL_STATIC.
//----------------------------------------------------------------------
#ifdef ARAP_STATIC
#define ARAP_DLL_API  // since ARAP_STATIC is defined, code is statically 
                      // linked and no exports or imports are needed
#else

#ifdef ARAP_DLL_EXPORTS
#define ARAP_DLL_API __declspec(dllexport)
#else
#define ARAP_DLL_API __declspec(dllimport)
#endif

#endif

#else
//----------------------------------------------------------------------
// ARAP_DLL_API is ignored for all other systems
//----------------------------------------------------------------------
#define ARAP_DLL_API
#endif

#endif