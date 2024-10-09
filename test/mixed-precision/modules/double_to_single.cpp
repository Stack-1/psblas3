#include <functional> // needed for bind


extern "C" void double_to_single(double *vect_double, float *vect_single, int size);


void double_to_single(float *vect_double, double *vect_single, int size){
    int i;

    for(i=0;i < size; i++){
        vect_single[i] = static_cast<float>(vect_double[i]);
    }
    
}