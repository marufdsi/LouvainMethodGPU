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

#ifndef COMMONCONSTANTS_H
#define	COMMONCONSTANTS_H

#include"devconstants.h"

#define CHUNK_PER_WARP 32
#define NR_THREAD_PER_BLOCK 128

#define WEIGHTED   0
#define UNWEIGHTED 1

#define PRINTALL 0

#endif	/* COMMONCONSTANTS_H */

