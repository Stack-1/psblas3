#include "cudadebug.h"
#include "cudalang.h"

extern "C"
{
#include "core.h"
#include "hell.h"
}

#include "debug.h"

#define VALUE_TYPE float
#define VALUE_TYPE_2 double
#define TYPE_SYMBOL S
#define TEX_FETCH_TYPE float
#include "hell_spmv_base_mx.cuh"