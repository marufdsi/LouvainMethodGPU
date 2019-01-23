/*

    Copyright (C) 2016, University of Bergen

    This file is part of Rundemanen - CUDA C++ parallel program for
    community detection

    Rundemanen is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rundemanen is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rundemanen.  If not, see <http://www.gnu.org/licenses/>.
    
    */

#ifndef OPENADDRESSING_H
#define	OPENADDRESSING_H

#include "hashitem.h"
#include"devconstants.h"

#ifdef RUNONGPU

__device__
#endif
inline unsigned int H1GPU(unsigned int key, unsigned int bucketSize) {
    return key % bucketSize;
}

#ifdef RUNONGPU

__device__
#endif
inline unsigned int H2GPU(unsigned int key, unsigned int bucketSize) {
    return 1 + (key % (bucketSize - 1));
}

#ifdef RUNONGPU
__device__
#endif
int hashInsertGPU(HashItem* Table, unsigned int* totNrAttempt, unsigned int bucketSize, HashItem *dataItem, float* tot, float wDegNode, float m2, float* bestGain, int* bestDest);

#ifdef RUNONGPU
__device__
#endif
int hashSearchGPU(HashItem* Table, unsigned int* totNrAttempt, unsigned int bucketSize, HashItem *dataItem);

#endif	/* OPENADDRESSING_H */

