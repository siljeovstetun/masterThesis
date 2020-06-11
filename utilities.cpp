#include "H5Cpp.h"
#include <iostream>
using namespace H5;


void writeToFileHDF5(double *x, double *y, double *angles, const int NP, const int NI, int L, double dt, int step, double *R){
  const H5std_string	FILE_NAME("/home/silje/Documents/masteroppgaveData/smitte/data.h5");
  const H5std_string	DATASET_NAME1("x");
  const H5std_string	DATASET_NAME2("y");
  const H5std_string	DATASET_NAME3("condition");
  const H5std_string  DATASET_NAME5("r");
  const H5std_string  DATASET_NAME4("info");

  const int	 NX = NP * (int)((NI + 1)/step);
  const int	 RANK = 1;

  double info[5] = {NP, NI, dt, L, step};
  H5File file(FILE_NAME, H5F_ACC_TRUNC);
	hsize_t dims[2];
	dims[0] = NX;
  dims[1] = 1;



	DataSpace dataspace(RANK, dims);
	DataSet dataset1 = file.createDataSet(DATASET_NAME1, PredType::NATIVE_DOUBLE, dataspace);
  DataSet dataset2 = file.createDataSet(DATASET_NAME2, PredType::NATIVE_DOUBLE, dataspace);
  DataSet dataset3 = file.createDataSet(DATASET_NAME3, PredType::NATIVE_DOUBLE, dataspace);

  dims[0] = (int)((NI + 1)/step);
  DataSpace dataspace2(RANK, dims);
  DataSet dataset5 = file.createDataSet(DATASET_NAME5, PredType::NATIVE_DOUBLE, dataspace2);


  dims[0] = 5;

  DataSpace dataspace1(RANK, dims);
  DataSet dataset4 = file.createDataSet(DATASET_NAME4, PredType::NATIVE_DOUBLE, dataspace1);
  dataset1.write(x, PredType::NATIVE_DOUBLE);
  dataset2.write(y, PredType::NATIVE_DOUBLE);
  dataset3.write(angles, PredType::NATIVE_DOUBLE);
  dataset4.write(info, PredType::NATIVE_DOUBLE);
  dataset5.write(R, PredType::NATIVE_DOUBLE);

}
