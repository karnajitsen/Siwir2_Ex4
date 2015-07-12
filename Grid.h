#include<iostream>
#include <assert.h>
#include <cmath>
#include <stdlib.h>
#include <memory>
#include <string.h>
#include<malloc.h>
#include "Cell.h"

class Grid{

	Cell * data;
	size_t xsize;
	size_t ysize;
	size_t zsize;
public:

	 Grid(size_t x, size_t y, size_t z)
	{
		xsize = x;
		ysize = y;
		zsize = z;
		data = new Cell[x * y * z];

		for (int i = 0; i < y; i++)
		{
			for (int j = 0; j < x; j++)
			{
				data[i*x + j].setZBoundary('D');
				data[i*x + j + x*y*(z-1)].setZBoundary('U');
			}
		}

		for (int i = 0; i < z; i++)
		{
			for (int j = 0; j < y; j++)
			{
				data[i*x*y + j*x].setXBoundary('W');
				data[i*x*y + j*x + x-1].setXBoundary('E');
			}
		}

		for (int i = 0; i < z; i++)
		{
			for (int j = 0; j < x; j++)
			{
				data[i*x*y + j].setYBoundary('S');
				data[i*x*y + (y - 1)*x +j].setYBoundary('N');
			}
		}
	}

	~Grid()
	{
		//--data;
		free(data);
	}

	inline Cell& operator()(const size_t x, const size_t y, const size_t z)
	{
		assert(x < xsize);
		assert(y < ysize);
		return data[x*xsize*ysize + y*xsize + z];
	}

	inline Cell& operator()(const size_t x, const size_t y, const size_t z) const
	{
		assert(x < xsize);
		assert(y < ysize);
		return data[x*xsize*ysize + y*xsize + z];
	}
};