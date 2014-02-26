#ifndef CUDAUTILITIES_H
#define CUDAUTILITIES_H

#define cudaErrorCheck(call) { cudaAssert(call,__FILE__,__LINE__); }
#define cudaFREE(ptr) if(ptr!=NULL) { cudaFree(ptr);ptr=NULL;}
void cudaAssert(const cudaError err, const char *file, const int line){}

#endif