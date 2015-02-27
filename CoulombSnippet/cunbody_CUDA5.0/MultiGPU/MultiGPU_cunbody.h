
#ifndef MULTIGPU_CUNBODY_H_INCLUDED
#define MULTIGPU_CUNBODY_H_INCLUDED

extern "C" {
    void cunbody1_force(double xj[][3], // position of j-th particles
                        double mj[],    // mass of j-th particles
                        double xi[][3], // position of i-th particles
                        double eps2,    // softening parameter
                        double ai[][3], // force of i-th particles
                        int ni,         // number of i-th particles
                        int nj);        // number of j-th particles
    void cunbody1_force_end(double xj[][3], // position of j-th particles
                            double mj[],    // mass of j-th particles
                            double xi[][3], // position of i-th particles
                            double eps2,    // softening parameter
                            double ai[][3], // force of i-th particles
                            int ni,         // number of i-th particles
                            int nj);        // number of j-th particles
    void cunbody1_force_begin(double xj[][3], // position of j-th particles
                              double mj[],    // mass of j-th particles
                              double xi[][3], // position of i-th particles
                              double eps2,    // softening parameter
                              double ai[][3], // force of i-th particles
                              int ni,         // number of i-th particles
                              int nj);        // number of j-th particles
}





void MultiGPUs_CUNBODY_begin();
void MultiGPUs_CUNBODY_end();
void CompaCPU();
#endif //  MULTIGPU_CUNBODY_H_INCLUDED
