#include <stdio.h>
#include <stdlib.h>


/* Remove elements by deleting from the lowest indices upward */
void remove_dead_cells_left(int *cellsX, int *numCells, int *cellsToKill, int numKillCells) {
    if (*numCells <= 0 || numKillCells <= 0) return;


    for (int i = 0; i < numKillCells; ++i) {
        /* Adjust for prior deletions: each previous delete shifts the array left by 1 */
        int dead = cellsToKill[i] - i;

        /* Shift left starting at the (adjusted) dead index */
        for (int j = dead; j < *numCells - 1; ++j) {
            cellsX[j] = cellsX[j + 1];
        }

        /* Optional: clear last slot */
        cellsX[*numCells - 1] = 0;

        /* Decrease logical size */
        (*numCells)--;
    }
}

int main(void) {
    /* Sanity test */
    int cellsX[10] = {6, 3, 5, 23, 4};
    int numCells = 5;

    int cellsToKill[2] = {0, 2};
    int numKillCells = 2;

    printf("Before: ");
    for (int i = 0; i < numCells; ++i) printf("%d ", cellsX[i]);
    printf("\n");

    remove_dead_cells_left(cellsX, &numCells, cellsToKill, numKillCells);

    printf("After:  ");
    for (int i = 0; i < numCells; ++i) printf("%d ", cellsX[i]);
    printf("\n");
    /* Expected: 3 23 4 */

    return 0;
}

