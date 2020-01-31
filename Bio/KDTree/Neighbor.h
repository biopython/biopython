/* Copyright 2008 by Michiel de Hoon.  All rights reserved.
 *
 * This file is part of the Biopython distribution and governed by your
 * choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
 * Please see the LICENSE file that should have been included as part of this
 * package.
 */

struct Neighbor;

struct Neighbor
{
        long int index1;
        long int index2;
        float radius;
	struct Neighbor* next;
};
