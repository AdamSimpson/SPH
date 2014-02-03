#ifndef fluid_renderer_h
#define fluid_renderer_h

void start_renderer();

// enum of displayed parameter values
typedef enum { 
    GRAVITY,
    VISCOSITY,
    DENSITY,
    PRESSURE,
    ELASTICITY   
 } parameters;

#endif
