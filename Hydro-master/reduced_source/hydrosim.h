#ifndef _HYDROSIM_H_
#define _HYDROSIM_H_

#include "LinkList.h"

void svSimulation(double dt,LinkList<2> &linklist);
void BSQSimulation(double dt,LinkList<2> &linklist);

template <int D>
void shear(LinkList<D>  &linklist);

template <int D>
void BSQshear(LinkList<D>  &linklist);


#endif
