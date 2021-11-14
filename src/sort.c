#include <stdio.h>
#include <stdlib.h>

/*INSERT SORT*/
void insertSort(int arr[], int arr2[], int arr3[], int r){
  int i,j,temp,temp2,temp3; 
	for(i=1;i<r;i++){ 
    temp=arr[i]; 
    temp2=arr2[i]; 
    temp3=arr3[i]; 
    j=i-1; 
    while(j>=0 && arr[j]>temp){
		  arr[j+1]=arr[j];
		  arr2[j+1]=arr2[j];
		  arr3[j+1]=arr3[j];
		  j--;		
	  }
	  arr[j+1]=temp;
	  arr2[j+1]=temp2;
	  arr3[j+1]=temp3;
	}
}

/* MERGE SORT */
/* https://www.geeksforgeeks.org/merge-sort/  
 * Merges two subarrays of arr[].
 * First subarray is arr[l..m], second subarray is arr[m+1..r] */
void mergeDouble(double arr[], int arr2[], int l, int m, int r){
  int i, j, k;
  int n1 = m - l + 1;
  int n2 = r - m;
  double *L,*R;
  int *L2,*R2;

  /* create temp arrays */
  L = malloc(n1*sizeof(double *));
  L2 = malloc(n1*sizeof(int *));
  R = malloc(n2*sizeof(double *));
  R2 = malloc(n2*sizeof(int *));

  /* Copy data to temp arrays L[] and R[] */
  for (i = 0; i < n1; i++){
    *(L+i) = arr[l + i];
    *(L2+i) = arr2[l + i];
  }
  for (j = 0; j < n2; j++){
    *(R+j) = arr[m + 1 + j];
    *(R2+j) = arr2[m + 1 + j];
  }

  /* Merge the temp arrays back into arr[l..r] */
  i = 0; // Initial index of first subarray
  j = 0; // Initial index of second subarray
  k = l; // Initial index of merged subarray
  while (i < n1 && j < n2){
    if (*(L+i) <= *(R+j)){
      arr[k] = *(L+i);
      arr2[k] = *(L2+i);
      i++;
    }
    else{
      arr[k] = *(R+j);
      arr2[k] = *(R2+j);
      j++;
      }
    k++;
  }

  /* Copy the remaining elements of L[], if there 
     are any */
  while (i < n1){
    arr[k] = *(L+i);
    arr2[k] = *(L2+i);
    i++;
    k++;
  }

  /* Copy the remaining elements of R[], if there 
     are any */
  while (j < n2){
    arr[k] = *(R+j);
    arr2[k] = *(R2+j);
    j++;
    k++;
  }
  free(L);
  free(L2);
  free(R);
  free(R2);
}


/* l is for left index and r is right index of the 
 * sub-array of arr to be sorted 
 * r == arr_size = sizeof(arr) / sizeof(arr[0]);
 * l == 0 */
void mergeSortDouble(double arr[],int arr2[], int l, int r){
  if (l < r) {
    /* Same as (l+r)/2, but avoids overflow for 
    * large l and h */
    int m = l + (r - l) / 2;

    /* Sort first and second halves */
    mergeSortDouble(arr,arr2, l, m);
    mergeSortDouble(arr,arr2, m + 1, r);
    mergeDouble(arr,arr2,l, m, r);
  }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void mergeFloat(float arr[], int arr2[], int l, int m, int r){
  int i, j, k;
  int n1 = m - l + 1;
  int n2 = r - m;
  float *L,*R;
  int *L2,*R2;

  /* create temp arrays */
  L = malloc(n1*sizeof(float *));
  L2 = malloc(n1*sizeof(int *));
  R = malloc(n2*sizeof(float *));
  R2 = malloc(n2*sizeof(int *));

  /* Copy data to temp arrays L[] and R[] */
  for (i = 0; i < n1; i++){
    *(L+i) = arr[l + i];
    *(L2+i) = arr2[l + i];
  }
  for (j = 0; j < n2; j++){
    *(R+j) = arr[m + 1 + j];
    *(R2+j) = arr2[m + 1 + j];
  }

  /* Merge the temp arrays back into arr[l..r] */
  i = 0; // Initial index of first subarray
  j = 0; // Initial index of second subarray
  k = l; // Initial index of merged subarray
  while (i < n1 && j < n2){
    if (*(L+i) <= *(R+j)){
      arr[k] = *(L+i);
      arr2[k] = *(L2+i);
      i++;
    }
    else{
      arr[k] = *(R+j);
      arr2[k] = *(R2+j);
      j++;
      }
    k++;
  }

  /* Copy the remaining elements of L[], if there 
     are any */
  while (i < n1){
    arr[k] = *(L+i);
    arr2[k] = *(L2+i);
    i++;
    k++;
  }

  /* Copy the remaining elements of R[], if there 
     are any */
  while (j < n2){
    arr[k] = *(R+j);
    arr2[k] = *(R2+j);
    j++;
    k++;
  }
  free(L);
  free(L2);
  free(R);
  free(R2);
}


/* l is for left index and r is right index of the 
 * sub-array of arr to be sorted 
 * r == arr_size = sizeof(arr) / sizeof(arr[0]);
 * l == 0 */
void mergeSortFloat(float arr[],int arr2[], int l, int r){
  if (l < r) {
    /* Same as (l+r)/2, but avoids overflow for 
    * large l and h */
    int m = l + (r - l) / 2;

    /* Sort first and second halves */
    mergeSortFloat(arr,arr2, l, m);
    mergeSortFloat(arr,arr2, m + 1, r);
    mergeFloat(arr,arr2,l, m, r);
  }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/* https://www.geeksforgeeks.org/merge-sort/ 
 * Merges two subarrays of arr[].
 * First subarray is arr[l..m], second subarray is arr[m+1..r] */
void mergeInt(int arr[], int arr2[], int l, int m, int r){
  int i, j, k;
  int n1 = m - l + 1;
  int n2 = r - m;
  int *L,*R;
  int *L2,*R2;

  /* create temp arrays */
  L = malloc(n1*sizeof(int *));
  L2 = malloc(n1*sizeof(int *));
  R = malloc(n2*sizeof(int *));
  R2 = malloc(n2*sizeof(int *));

  /* Copy data to temp arrays L[] and R[] */
  for (i = 0; i < n1; i++){
    *(L+i) = arr[l + i];
    *(L2+i) = arr2[l + i];
  }
  for (j = 0; j < n2; j++){
    *(R+j) = arr[m + 1 + j];
    *(R2+j) = arr2[m + 1 + j];
  }

  /* Merge the temp arrays back into arr[l..r] */
  i = 0; // Initial index of first subarray
  j = 0; // Initial index of second subarray
  k = l; // Initial index of merged subarray
  while (i < n1 && j < n2){
    if (*(L+i) <= *(R+j)){
      arr[k] = *(L+i);
      arr2[k] = *(L2+i);
      i++;
    }
    else{
      arr[k] = *(R+j);
      arr2[k] = *(R2+j);
      j++;
      }
    k++;
  }

  /* Copy the remaining elements of L[], if there 
     are any */
  while (i < n1){
    arr[k] = *(L+i);
    arr2[k] = *(L2+i);
    i++;
    k++;
  }

  /* Copy the remaining elements of R[], if there 
     are any */
  while (j < n2){
    arr[k] = *(R+j);
    arr2[k] = *(R2+j);
    j++;
    k++;
  }
  free(L);
  free(L2);
  free(R);
  free(R2);
}


/* l is for left index and r is right index of the
 * sub-array of arr to be sorted
 * r == arr_size = sizeof(arr) / sizeof(arr[0]);
 * l == 0 */
void mergeSortInt(int arr[],int arr2[], int l, int r){
  if (l < r) {
    /* Same as (l+r)/2, but avoids overflow for
    * large l and h */
    int m = l + (r - l) / 2;

    /* Sort first and second halves */
    mergeSortInt(arr,arr2, l, m);
    mergeSortInt(arr,arr2, m + 1, r);
    mergeInt(arr,arr2,l, m, r);
  }
}



/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/* https://www.geeksforgeeks.org/merge-sort/ 
 * Merges two subarrays of arr[].
 * First subarray is arr[l..m], second subarray is arr[m+1..r] */
void mergeInt3(int arr[], int arr2[], int arr3[], int l, int m, int r){
  int i, j, k;
  int n1 = m - l + 1;
  int n2 = r - m;
  int *L,*R;
  int *L2,*R2;
  int *L3,*R3;

  /* create temp arrays */
  L = malloc(n1*sizeof(int *));
  L2 = malloc(n1*sizeof(int *));
  L3 = malloc(n1*sizeof(int *));
  R = malloc(n2*sizeof(int *));
  R2 = malloc(n2*sizeof(int *));
  R3 = malloc(n2*sizeof(int *));

  /* Copy data to temp arrays L[] and R[] */
  for (i = 0; i < n1; i++){
    *(L+i) = arr[l + i];
    *(L2+i) = arr2[l + i];
    *(L3+i) = arr3[l + i];
  }
  for (j = 0; j < n2; j++){
    *(R+j) = arr[m + 1 + j];
    *(R2+j) = arr2[m + 1 + j];
    *(R3+j) = arr3[m + 1 + j];
  }

  /* Merge the temp arrays back into arr[l..r] */
  i = 0; // Initial index of first subarray
  j = 0; // Initial index of second subarray
  k = l; // Initial index of merged subarray
  while (i < n1 && j < n2){
    if (*(L+i) <= *(R+j)){
      arr[k] = *(L+i);
      arr2[k] = *(L2+i);
      arr3[k] = *(L3+i);
      i++;
    }
    else{
      arr[k] = *(R+j);
      arr2[k] = *(R2+j);
      arr3[k] = *(R3+j);
      j++;
      }
    k++;
  }

  /* Copy the remaining elements of L[], if there 
     are any */
  while (i < n1){
    arr[k] = *(L+i);
    arr2[k] = *(L2+i);
    arr3[k] = *(L3+i);
    i++;
    k++;
  }

  /* Copy the remaining elements of R[], if there 
     are any */
  while (j < n2){
    arr[k] = *(R+j);
    arr2[k] = *(R2+j);
    arr3[k] = *(R3+j);
    j++;
    k++;
  }
  free(L);
  free(L2);
  free(L3);
  free(R);
  free(R2);
  free(R3);
}


/* l is for left index and r is right index of the
 * sub-array of arr to be sorted
 * r == arr_size = sizeof(arr) / sizeof(arr[0]);
 * l == 0 */
void mergeSortInt3(int arr[],int arr2[], int arr3[], int l, int r){
  if (l < r) {
    /* Same as (l+r)/2, but avoids overflow for
    * large l and h */
    int m = l + (r - l) / 2;

    /* Sort first and second halves */
    mergeSortInt3(arr, arr2, arr3, l, m);
    mergeSortInt3(arr, arr2, arr3, m + 1, r);
    mergeInt3(arr,arr2, arr3,l, m, r);
  }
}
