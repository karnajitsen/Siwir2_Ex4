#include<iostream>
#include<fstream>
using namespace std;

inline void writeVTK(string filename, Grid *grd)
{
	//ofstream file;
	////cout << filename;
	//file.open(filename);
	//int len = (grd->getXsize() - 2) * (grd->getYsize() - 2);
	//file << "# vtk DataFile Version 4.0" << std::endl;
	//file << "SiwiRVisFile" << std::endl;
	//file << "ASCII" << std::endl;
	//file << "DATASET STRUCTURED_POINTS" << std::endl;
	//file << "DIMENSIONS " << grd->getXsize()-2 << " " << grd->getYsize()-2 << " 1"  << std::endl;
	//file << "ORIGIN 0 0 0" << std::endl;
	//file << "POINT_DATA " << len << std::endl;
	//file << std::endl;
	//file << std::endl;
	//file << "SCALARS flags double 1" << len << std::endl;
	//file << "LOOKUP_TABLE default" << len << std::endl;
	//
	//for (int i = 0; i < len; i++)
	//{
	//	file << "1" << std::endl;
	//}
	//file << std::endl;

	//file << "SCALARS density double 1" << std::endl;

	//file << "LOOKUP_TABLE default" << std::endl;

	//for (size_t i = 1; i < grd->getYsize()-1; i++)
	//{
	//	for (size_t j = 1; j < grd->getXsize()-1; j++)
	//	{
	//		file << (*grd)(i,j,11) << std::endl;
	//	}
	//}

	//file << std::endl;

	//file << "VECTORS velocity double" << std::endl;

	//for (size_t i = 1; i < grd->getYsize() - 1; i++)
	//{
	//	for (size_t j = 1; j < grd->getXsize() - 1; j++)
	//	{
	//		file << (*grd)(i, j, 9) << " " << (*grd)(i, j, 10) << " 0" << std::endl;
	//	}
	//}
	//file << std::endl;
	//file << std::endl;
	//file.close();
}
