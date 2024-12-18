// point to IP
const int nWireGroup_1 = 14;
double w_rows_1[nWireGroup_1] = {0,  0,0, 1,1,2,3,3,4,4,4,5,5,5};
double w_cols_1[nWireGroup_1] = {-2,-1,0,-1,0,0,0,1,0,1,2,0,1,2};
double w_data_1[nWireGroup_1] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1};

const int nWireGroup_2 = 13;
double w_rows_2[nWireGroup_2] = {0,0,0,1,1,2,3,3,4,4,5,5,5};
double w_cols_2[nWireGroup_2] = {-3,-2,-1,-1,0,0,0,1,1,2,1,2,3};
double w_data_2[nWireGroup_2] = {1,1,1,1,1,1,1,1,1,1,1,1,1};

const int nDefinedWirePatterns = 2; 

int nWGsInPatterns[nDefinedWirePatterns] = {nWireGroup_1, nWireGroup_2};
int patternRanks_w[nDefinedWirePatterns] = {1, 2};

double* w_rows[nDefinedWirePatterns] = {w_rows_1, w_rows_2};
double* w_cols[nDefinedWirePatterns] = {w_cols_1, w_cols_2};
double* w_data[nDefinedWirePatterns] = {w_data_1, w_data_2};


/*
// Keep only following ones
-----------111--------------------------------------------------------------------------------------------------
-----------11---------------------------------------------------------------------------------------------------
-----------1----------------------------------------------------------------------------------------------------
----------11----------------------------------------------------------------------------------------------------
---------111----------------------------------------------------------------------------------------------------
----------11----------------------------------------------------------------------------------------------------


------------111-------------------------------------------------------------------------------------------------
-----------11---------------------------------------------------------------------------------------------------
-----------1----------------------------------------------------------------------------------------------------
----------11----------------------------------------------------------------------------------------------------
---------11-----------------------------------------------------------------------------------------------------
--------111-----------------------------------------------------------------------------------------------------
*/
