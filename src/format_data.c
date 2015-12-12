/*
 *
 * format_data.c 
 *
 * -----------------------------------------------------------------------------
 *
 * Global variable definitions:
 *      None.
 *
 * -----------------------------------------------------------------------------
 *
 * Subroutines: 
 *      **read_data_file
 *
 */

#include "format_data.h"
#include "deconvolution_main.h"

/*
 * read_data_file: this scans the inputted data file and saves the second two
 *                 columns; they contain the time and the concentration. It
 *                 log-transforms the concentration before outputting the data
 *
 *      ARGUMENTS: 
 *          char *datafile; the name of inputted data file
 *          int *N; the variable where we will store the number of observations
 *
 *      RETURNS: 
 *          data; the Nx2 matrix of times and log-concentrations
 *
 */
/*{{{*/
/*
 * -----------------------------------------------------------------------------
 *
 * Variable definitions
 *      i       - generic counter
 *      j       - where we store the observation number
 *      *labels - this is where we store the column names
 *      x       - we store the data in column 2 here while we count the number of rows
 *      y       - we store the data in column 3 here while we count the number of rows
 *      **data  - we read the entries in columns 2 and 3 (minus the header) into the matrix data
 *      *inf    - infile variable
 *
 * -----------------------------------------------------------------------------
 *
 * Subroutines used
 *      None
 *
 */

double **read_data_file(char *datafile, int *N) {

  int i;           // Generic counter
  int j;           // Stores observation number
  int placeholder; // Placeholder for fscanf's return value (not used otherwise)
  char *labels;    // Column names
  double x;        // holder for column 2 (time)
  double y;        // holder for column 3 (concentration)
  double **data;   // Nx2 matrix of time and log-concentration
  FILE *inf;       // Infile variable

  inf = fopen(datafile, "r");

  if (inf == NULL) {
    Rprintf("file %s does not exist\n", datafile);
    free(N);
    exit(0);
  }

  *N = 0;

  // not interested in the column names, read them in and forget about them 
  labels = (char *)calloc(50, sizeof(char *));
  placeholder = fscanf(inf, "%s %s %s\n", labels, labels, labels);

  // Count number of lines
  while (fscanf(inf, "%d %lf %lf\n", &j, &x, &y) != EOF) {
    (*N)++;
  }

  // Move back to start of file and store/remove labels
  rewind(inf);
  placeholder = fscanf(inf, "%s %s %s\n", labels, labels, labels);
  free(labels);

  //--------------------------------------------------------
  // Allocate memory to hold the hormonal time series and read into an Nx2
  // matrix.  
  //
  // each row contains a recording of the time series.  The first column is
  // the time the measurement was taken, the second column records the log
  // hormonal concentration (note: we are not using the obs. #) 
  //--------------------------------------------------------
  data = (double **)calloc(*N, sizeof(double *));
  for (i = 0; i < *N; i++) {
    data[i] = (double *)calloc(2, sizeof(double));
  }

  // Read the time series into memory                          
  for (i = 0; i < *N; i++) {
    placeholder = fscanf(inf, "%d %lf %lf\n", &j, &data[i][0], &data[i][1]);
    data[i][1] = log(data[i][1]);
  }

  fclose(inf);

  // Return data matrix 
  return data;

}
/*}}}*/




/*******************************************************************************
 * END OF FILE
 ******************************************************************************/

