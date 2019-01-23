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

#ifndef HOSTCONSTANTS_H
#define	HOSTCONSTANTS_H


#define SHARED_TABLE_SIZE 479 // Must be a Prime
#define WARP_TABLE_SIZE_1 127 // Must be a Prime

#define CAPACITY_FACTOR_NUMERATOR 2
#define CAPACITY_FACTOR_DENOMINATOR 3
#define HALF_WARP 16
#define QUARTER_WARP 8
#endif	/* HOSTCONSTANTS_H */
