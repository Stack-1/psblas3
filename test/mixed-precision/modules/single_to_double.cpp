#include <functional> // needed for bind

extern "C" void single_to_double(float *vect_single, double *vect_double, int size);


void single_to_double(float *vect_single, double *vect_double, int size){
    int i;

    for(i=0;i < size; i++){
        vect_double[i] = static_cast<double>(vect_single[i]);
    }
    
}