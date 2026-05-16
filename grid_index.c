#include "grid_index.h"

int grid_index(int N_S, int ix, int iy)
/*
   Goal
   ===
   This function calculates, in lexicographical order, the index of the equation according to the position in the grid.
   
   Arguments
   ===========
   (input)   N_S   - Number of point in the grid in S (x-axis)
   (input)   ix    - Current value of x in the grid
   (input)   iy    - Current value of y in the grid
                
   (output) ind    - Associated index
*/
{
    int ind = ix + N_S*iy;

    return ind;
}