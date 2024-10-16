const int nHalfStrip_1 = 12;
double s_rows_1[nHalfStrip_1] = {0,0,0,1,2,3,4,4,4,5,5,5};
double s_cols_1[nHalfStrip_1] = {-1,0,1,0,0,0,-1,0,1,-1,0,1};
double s_data_1[nHalfStrip_1] = {1,1,1,1,1,1,1,1,1,1,1,1};

const int nHalfStrip_2 = 14;
double s_rows_2[nHalfStrip_2] = {0,0,0,1,1,2,3,3,4,4,4,5,5,5};
double s_cols_2[nHalfStrip_2] = {0,1,2,0,1,0,-1,0,-2,-1,0,-2,-1,0};
double s_data_2[nHalfStrip_2] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1};

const int nHalfStrip_3 = 14;
double s_rows_3[nHalfStrip_3] = {0,0,0,1,1,2,3,3,4,4,4,5,5,5};
double s_cols_3[nHalfStrip_3] = {-2,-1,0,-1,0,0,0,1,0,1,2,0,1,2};
double s_data_3[nHalfStrip_3] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1};

const int nHalfStrip_4 = 13;
double s_rows_4[nHalfStrip_4] = {0,0,0,1,1,2,3,3,4,4,5,5,5};
double s_cols_4[nHalfStrip_4] = {1,2,3,0,1,0,-1,0,-2,-1,-3,-2,-1};
double s_data_4[nHalfStrip_4] = {1,1,1,1,1,1,1,1,1,1,1,1,1};

const int nHalfStrip_5 = 13;
double s_rows_5[nHalfStrip_5] = {0,0,0,1,1,2,3,3,4,4,5,5,5};
double s_cols_5[nHalfStrip_5] = {-3,-2,-1,-1,0,0,0,1,1,2,1,2,3};
double s_data_5[nHalfStrip_5] = {1,1,1,1,1,1,1,1,1,1,1,1,1};

const int nHalfStrip_6 = 15;
double s_rows_6[nHalfStrip_6] = {0,0,0,0,1,1,2,3,3,4,4,4,5,5,5};
double s_cols_6[nHalfStrip_6] = {1,2,3,4,1,2,0,-2,-1,-4,-3,-2,-4,-3,-2};
double s_data_6[nHalfStrip_6] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

const int nHalfStrip_7 = 15;
double s_rows_7[nHalfStrip_7] = {0,0,0,0,1,1,2,3,3,4,4,4,5,5,5};
double s_cols_7[nHalfStrip_7] = {-4,-3,-2,-1,-2,-1,0,1,2,2,3,4,2,3,4};
double s_data_7[nHalfStrip_7] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

const int nHalfStrip_8 = 15;
double s_rows_8[nHalfStrip_8] = {0,0,0,1,1,2,3,3,3,4,4,4,5,5,5};
double s_cols_8[nHalfStrip_8] = {3,4,5,1,2,0,-2,-1,0,-4,-3,-2,-5,-4,-3};
double s_data_8[nHalfStrip_8] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

const int nHalfStrip_9 = 15;
double s_rows_9[nHalfStrip_9] = {0,0,0,1,1,2,3,3,3,4,4,4,5,5,5};
double s_cols_9[nHalfStrip_9] = {-5,-4,-3,-2,-1,0,0,1,2,2,3,4,3,4,5};
double s_data_9[nHalfStrip_9] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

// above are from NIMA_592(2008)26

const int nHalfStrip_10 = 18;
double s_rows_10[nHalfStrip_10] = {0,0,0,0,1,1,1,2,3,3,3,4,4,4,5,5,5,5};
double s_cols_10[nHalfStrip_10] = {-6,-5,-4,-3,-3,-2,-1,0,1,2,3,4,5,6,6,7,8,9};
double s_data_10[nHalfStrip_10] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

const int nHalfStrip_11 = 18;
double s_rows_11[nHalfStrip_11] = {0,0,0,0,1,1,1,2,3,3,3,4,4,4,5,5,5,5};
double s_cols_11[nHalfStrip_11] = {3,4,5,6,1,2,3,0,-3,-2,-1,-6,-5,-4,-9,-8,-7,-6};
double s_data_11[nHalfStrip_11] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

int nHalfStrips[9] = {12,14,14,13,13,15,15,15,15};
int patternRanks[9] = {1,2,2,3,3,4,4,5,5};
double* s_rows[9] = {s_rows_1, s_rows_2, s_rows_3, s_rows_4, s_rows_5, s_rows_6, s_rows_7, s_rows_8, s_rows_9};
double* s_cols[9] = {s_cols_1, s_cols_2, s_cols_3, s_cols_4, s_cols_5, s_cols_6, s_cols_7, s_cols_8, s_cols_9};
double* s_data[9] = {s_data_1, s_data_2, s_data_3, s_data_4, s_data_5, s_data_6, s_data_7, s_data_8, s_data_9};
