#include <stdio.h>
#include <stdlib.h>

/* https://www.geeksforgeeks.org/binary-search/ */

//recursive
/*int binarySearch( int arr[], int l, int r, int id){
  if (r >= l) {
    int mid = l + (r - l) / 2;

    // If the element is present at the middle
    // itself
    if (arr[mid] == id)
      return mid;

    // If element is smaller than mid, then
    // it can only be present in left subarray
    if (arr[mid] > id)
      return binarySearch(arr, l, mid - 1, id);

    // Else the element can only be present
    // in right subarray
    return binarySearch(arr, mid + 1, r, id);
    }

    // We reach here when element is not
    // present in array
    printf("halo tag of the particle is not found in the catalogue\n");
    exit(42);
}*/

//iterative
int binarySearch(int arr[], int l, int r, int x){
    while (l <= r) {
        int m = l + (r - l) / 2;

        // Check if x is present at mid
        if (arr[m] == x)
            return m;

        // If x greater, ignore left half
        if (arr[m] < x)
            l = m + 1;

        // If x is smaller, ignore right half
        else
            r = m - 1;
    }
    printf("halo tag of the particle is not found in the catalogue\n");
    exit(42);
    // if we reach here, then element was
    // not present
}
