#include <iostream>
#include<fstream>
#include<cmath>
#include<memory>
#include<stdlib.h>

#include "FileReader.h"
#include "VTKFileWriter.h"
#include "Grid.h"

#define DIM 3
typedef double Real;

size_t vis_space;
Real t_start;
Real t_end;
Real delta_t;
Real x_min;
Real y_min;
Real z_min;
Real x_max;
Real y_max;
Real z_max;
Real r_cut;
Real epsilon;
Real sigma;
size_t xsize;
size_t ysize;
size_t zsize;
size_t cellx;
size_t celly;
size_t cellz;
Grid * grd = nullptr;
Grid * tgrd = nullptr;

void getCellNumber(Real x, Real y, Real z, size_t& i, size_t& j, size_t& k)
{
	i = (size_t)ceil(fmod(x , x_max) / cellx);
	j = (size_t)ceil(fmod(y, y_max) / celly);
	k = (size_t)ceil(fmod(z, z_max) / cellz);
}
void inline init()
{
	grd = new Grid( xsize , ysize , zsize);
	tgrd = new Grid(xsize, ysize, zsize);
	
	ifstream blocks;
	blocks.open("blocks-small.dat");
	Real a, b, c, d, e, f, g;
	size_t i, j, k;
	//i = 0, j = 0, k = 0;
	while (blocks >> a && blocks >> b && blocks >> c && blocks >> d && blocks >> e && blocks >> f && blocks >> g)
	{
		Particle p;
		p.m = a;
		p.x[0] = b;
		p.x[1] = c;
		p.x[2] = d;
		p.v[0] = e;
		p.v[1] = f;
		p.v[2] = g;
		p.F[0] = 0.0;
		p.F[1] = 0.0;
		p.F[2] = 0.0;
		getCellNumber(b, c, d, i, j, k);
		(*grd)(i, j, k).addParticle(p);

	}	
}

inline Real compRij2(Particle a, Particle b)
{
	return ((a.x[0] - b.x[0]) * (a.x[0] - b.x[0]) + (a.x[1] - b.x[1]) * (a.x[1] - b.x[1]) + (a.x[2] - b.x[2]) * (a.x[2] - b.x[2]));
}

inline void compForceCenterCell(size_t z, size_t y, size_t x, size_t p)
{
	size_t len = (*grd)(z, y, x).getLength();
	//Particle cp = (*grd)(z, y, x)(p);
	//Particle np;
	Real f,s;
	double rij2 = 0.0;
	for (size_t i = 0; i < len || i!=p; i++)
	{
		rij2 = compRij2((*grd)(z, y, x)(p), (*grd)(z, y, x)(i));
		s = sigma * sigma / rij2;
		s = s * s * s;
		f = 24.0 * epsilon * s * (1.0 - 2.0 * s) / rij2;
		for (size_t d = 0; d < DIM ; d++)
		(*grd)(z, y, x)(i).F[d] += f * ((*grd)(z, y, x)(i).x[d] - (*grd)(z, y, x)(p).x[d]);
	}
}

inline void compForceBetCells(size_t z, size_t y, size_t x, size_t nz, size_t ny, size_t nx,size_t p)
{
	size_t len = (*grd)(nz,ny, nx).getLength();
	
	Real f, s;
	double rij2 = 0.0;
	for (size_t i = 0; i < len || i != p; i++)
	{
		rij2 = compRij2((*grd)(z, y, x)(p), (*grd)(nz, ny, nx)(i));
		if (sqrt(rij2) > r_cut)
			continue;
		s = sigma * sigma / rij2;
		s = s * s * s;
		f = 24.0 * epsilon * s * (1.0 - 2.0 * s) / rij2;
		for (size_t d = 0; d < DIM; d++)
			(*grd)(z, y, x)(i).F[d] += f * ((*grd)(z, y, x)(i).x[d] - (*grd)(nz, ny, nx)(p).x[d]);
	}
}

inline void compForceUpCells(size_t z, size_t y, size_t x, size_t p)
{
	size_t len = (*grd)(z, y, x).getLength();
	size_t uz, ex, wx, sy, ny;
	if ((*grd)(z, y, x).getZBoundary() == 'U')
		uz = 0;
	else
		uz = z + 1;
	if ((*grd)(z, y, x).getXBoundary() == 'E')
		ex = 0;
	else
		ex = x + 1;

	if ((*grd)(z, y, x).getXBoundary() == 'W')
		wx = xsize - 1;
	else
		wx = x - 1;

	if ((*grd)(z, y, x).getYBoundary() == 'S')
		sy = ysize - 1;
	else
		sy = y - 1;

	if ((*grd)(z, y, x).getYBoundary() == 'N')
		ny = 0;
	else
		ny = y + 1;
	
	compForceBetCells(x,y,z,x,y,uz,p);
	compForceBetCells(x, y, z, ex, y, uz,p);
	compForceBetCells(x, y, z, wx, y, uz, p);
	compForceBetCells(x, y, z, x, sy, uz, p);
	compForceBetCells(x, y, z, x, ny, uz, p);
	compForceBetCells(x, y, z, ex, sy, uz, p);
	compForceBetCells(x, y, z, ex, ny, uz, p);
	compForceBetCells(x, y, z, wx, sy, uz, p);
	compForceBetCells(x, y, z, wx, ny, uz, p);
	
}


inline void compForceDownCells(size_t z, size_t y, size_t x, size_t p)
{
	size_t len = (*grd)(z, y, x).getLength();
	size_t uz, ex, wx, sy, ny;
	if ((*grd)(z, y, x).getZBoundary() == 'D')
		uz = zsize -1;
	else
		uz = z - 1;
	if ((*grd)(z, y, x).getXBoundary() == 'E')
		ex = 0;
	else
		ex = x + 1;

	if ((*grd)(z, y, x).getXBoundary() == 'W')
		wx = xsize - 1;
	else
		wx = x - 1;

	if ((*grd)(z, y, x).getYBoundary() == 'S')
		sy = ysize - 1;
	else
		sy = y - 1;

	if ((*grd)(z, y, x).getYBoundary() == 'N')
		ny = 0;
	else
		ny = y + 1;

	compForceBetCells(x, y, z, x, y, uz, p);
	compForceBetCells(x, y, z, ex, y, uz, p);
	compForceBetCells(x, y, z, wx, y, uz, p);
	compForceBetCells(x, y, z, x, sy, uz, p);
	compForceBetCells(x, y, z, x, ny, uz, p);
	compForceBetCells(x, y, z, ex, sy, uz, p);
	compForceBetCells(x, y, z, ex, ny, uz, p);
	compForceBetCells(x, y, z, wx, sy, uz, p);
	compForceBetCells(x, y, z, wx, ny, uz, p);

}


inline void compForceEastCells(size_t z, size_t y, size_t x, size_t p)
{
	size_t len = (*grd)(z, y, x).getLength();
	size_t  ex, sy, ny;
	
	if ((*grd)(z, y, x).getXBoundary() == 'E')
		ex = 0;
	else
		ex = x + 1;

	if ((*grd)(z, y, x).getYBoundary() == 'S')
		sy = ysize - 1;
	else
		sy = y - 1;

	if ((*grd)(z, y, x).getYBoundary() == 'N')
		ny = 0;
	else
		ny = y + 1;

	compForceBetCells(x, y, z, ex, y, z, p);
	compForceBetCells(x, y, z, ex, sy, z, p);
	compForceBetCells(x, y, z, ex, ny, z, p);
	
}


inline void compForceWestCells(size_t z, size_t y, size_t x, size_t p)
{
	size_t len = (*grd)(z, y, x).getLength();
	size_t  wx, sy, ny;
	
	if ((*grd)(z, y, x).getXBoundary() == 'W')
		wx = xsize -1;
	else
		wx = x - 1;

	if ((*grd)(z, y, x).getYBoundary() == 'S')
		sy = ysize - 1;
	else
		sy = y + 1;

	if ((*grd)(z, y, x).getYBoundary() == 'N')
		ny = 0;
	else
		ny = y - 1;

	compForceBetCells(x, y, z, wx, y, z, p);
	compForceBetCells(x, y, z, wx, sy, z, p);
	compForceBetCells(x, y, z, wx, ny, z, p);

}


inline void compForceSouthCells(size_t z, size_t y, size_t x, size_t p)
{
	size_t len = (*grd)(z, y, x).getLength();
	size_t  sy;

	if ((*grd)(z, y, x).getYBoundary() == 'S')
		sy = ysize - 1;
	else
		sy = y - 1;

	compForceBetCells(x, y, z, x, sy, z, p);
	
}

inline void compForceNorthCells(size_t z, size_t y, size_t x, size_t p)
{
	size_t len = (*grd)(z, y, x).getLength();
	size_t  ny;

	if ((*grd)(z, y, x).getYBoundary() == 'N')
		ny = 0;
	else
		ny = y + 1;

	compForceBetCells(x, y, z, x, ny, z, p);

}

inline void compForce()
{
	//Particle ** np;
	for (size_t i = 0; i < zsize; i++)
	{
		for (size_t j = 0; j < ysize; j++)
		{
			for (size_t k = 0; k < xsize; k++)
			{
				size_t len = (*grd)(i, j, k).getLength();
				//findRcutParticles(i, j, k, np);
				for (size_t p = 0; p < len; p++)
				{
					(*grd)(i, j, k)(p).FOld[0] = (*grd)(i, j, k)(p).F[0];
					(*grd)(i, j, k)(p).FOld[1] = (*grd)(i, j, k)(p).F[1];
					(*grd)(i, j, k)(p).FOld[2] = (*grd)(i, j, k)(p).F[2];

					compForceCenterCell(i, j, k,p);
					compForceUpCells(i, j, k,p);
					compForceDownCells(i, j, k, p);
					compForceSouthCells(i, j, k, p);
					compForceNorthCells(i, j, k, p);
					compForceEastCells(i, j, k, p);
					compForceWestCells(i, j, k, p);
				}
			}
		}
	}
}

inline void compPosition()
{
	
	size_t a, b, c;
	for (size_t i = 0; i < zsize; i++)
	{
		for (size_t j = 0; j < ysize; j++)
		{
			for (size_t k = 0; k < xsize; k++)
			{
				size_t len = (*grd)(i, j, k).getLength();
				for (size_t p = 0; p < len; p++)
				{
					(*grd)(i, j, k)(p).x[0] = (*grd)(i, j, k)(p).x[0] + delta_t * (*grd)(i, j, k)(p).v[0] + (*grd)(i, j, k)(p).F[0] * delta_t * delta_t * 0.5 / (*grd)(i, j, k)(p).m;
					(*grd)(i, j, k)(p).x[1] = (*grd)(i, j, k)(p).x[1] + delta_t * (*grd)(i, j, k)(p).v[1] + (*grd)(i, j, k)(p).F[1] * delta_t * delta_t * 0.5 / (*grd)(i, j, k)(p).m;
					(*grd)(i, j, k)(p).x[2] = (*grd)(i, j, k)(p).x[2] + delta_t * (*grd)(i, j, k)(p).v[2] + (*grd)(i, j, k)(p).F[2] * delta_t * delta_t * 0.5 / (*grd)(i, j, k)(p).m;
					getCellNumber((*grd)(i, j, k)(p).x[0], (*grd)(i, j, k)(p).x[1], (*grd)(i, j, k)(p).x[2], a, b, c);
					if (a != k || b != j || c != i)
					{
						(*grd)(c, b, a).addParticle((*grd)(i, j, k)(p));
						(*grd)(i, j, k).removeParticle(p);
					}
				}
			}
		}
	}
}

inline void compVelocity()
{
	for (size_t i = 0; i < zsize; i++)
	{
		for (size_t j = 0; j < ysize; j++)
		{
			for (size_t k = 0; k < xsize; k++)
			{
				size_t len = (*grd)(i, j, k).getLength();
				for (size_t p = 0; p < len; p++)
				{
					(*grd)(i, j, k)(p).v[0] += delta_t * 0.5 * ((*grd)(i, j, k)(p).F[0] + (*grd)(i, j, k)(p).FOld[0]) / (*grd)(i, j, k)(p).m;
					(*grd)(i, j, k)(p).v[1] += delta_t * 0.5 * ((*grd)(i, j, k)(p).F[1] + (*grd)(i, j, k)(p).FOld[1]) / (*grd)(i, j, k)(p).m;
					(*grd)(i, j, k)(p).v[2] += delta_t * 0.5 * ((*grd)(i, j, k)(p).F[2] + (*grd)(i, j, k)(p).FOld[2]) / (*grd)(i, j, k)(p).m;
				}
			}
		}
	}
}

inline void performTimeSteps()
{
	Real t = t_start;
	compForce();
	t = t + delta_t;
	while (t <= t_end)
	{
		compPosition();
		compForce();
		compVelocity();
		//compThermodynamics();
	}
}

int main(int argc, char** argv)
{

	if (argc < 2)
	{
		std::cout << "Invalid number of argument";
		exit(0);
	}

	string fname = argv[1];
	ifstream paramfile;
	string tmp;
	string vtkfilename;
	size_t vtk_step;

	FileReader* fr = new FileReader();

	fr->readParameters(fname);

	vis_space = fr->getParameter<size_t>("vis_space");
	t_start = fr->getParameter<Real>("t_start");
	t_end = fr->getParameter<Real>("t_end");
	delta_t = fr->getParameter<Real>("delta_t");
	x_min = fr->getParameter<Real>("x_min");
	y_min = fr->getParameter<Real>("y_min");
	z_min = fr->getParameter<Real>("z_min");
	x_max = fr->getParameter<Real>("x_max");
	y_max = fr->getParameter<Real>("y_max");
	z_max = fr->getParameter<Real>("z_max");
	r_cut = fr->getParameter<Real>("r_cut");
	epsilon = fr->getParameter<Real>("epsilon");
	sigma = fr->getParameter<Real>("sigma");
	cellx = r_cut;
	celly = r_cut;
	cellz = r_cut;
	xsize = (size_t)ceil((x_max - x_min) / r_cut);
	ysize = (size_t)ceil((y_max - y_min) / r_cut);
	zsize = (size_t)ceil((z_max - z_min) / r_cut);
	init();
	performTimeSteps();

}