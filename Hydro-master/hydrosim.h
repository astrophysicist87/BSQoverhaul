#ifndef _HYDROSIM_H_
#define _HYDROSIM_H_

#include "LinkList.h"

void Simulation(double dt,LinkList<2> &linklist);
void svSimulation(double dt,LinkList<2> &linklist);

template <int D> 
void idealhydro3(LinkList<D> &linklist);

template <int D> 
void shear(LinkList<D>  &linklist);

#endif
