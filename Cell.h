#include <iostream>
#include<fstream>
#include<cmath>
#include<memory>
#include<stdlib.h>
#define DIM 3
typedef double Real;

typedef struct {
	Real m;
	Real x[DIM];
	Real v[DIM];
	Real F[DIM];
	Real FOld[DIM];
}Particle;

typedef struct ParticleNode{
	Particle p;
	struct ParticleNode * next;
}ParticleNode;

class Cell {

	ParticleNode * header;
	size_t length;
	char xbndy, ybndy, zbndy;
	public:
	Cell()
	{
		header = nullptr;
		length = 0;
		xbndy = 'C';
		ybndy = 'C';
		zbndy = 'C';
	}

	inline void setXBoundary(char c)
	{
		xbndy = c;
	}

	inline void setYBoundary(char c)
	{
		ybndy = c;
	}

	inline void setZBoundary(char c)
	{
		zbndy = c;
	}

	inline char getXBoundary()
	{
		return xbndy;
	}

	inline char getYBoundary()
	{
		return ybndy;
	}

	inline char getZBoundary()
	{
		return zbndy;
	}

	inline size_t getLength()
	{
		return length;
	}
	inline Particle& operator()(size_t i, Particle cp)
	{
		assert(i < length);
		ParticleNode * p = header;
		int j;
		for (j = 0; j < length; j++)
		{
			if (j == i && (p->p.x[0] != cp.x[0] || p->p.x[1] != cp.x[1] || p->p.x[2] != cp.x[2]))
				return p->p;
			else if (j == i)
				return p->next->p;
			else
				p = p->next;
		}
	}

	inline Particle& operator()(size_t i, Particle cp) const
	{
		assert(i < length);
		ParticleNode * p = header;
		int j;
		for (j = 0; j < length; j++)
		{
			if (j == i && (p->p.x[0] != cp.x[0] || p->p.x[1] != cp.x[1] || p->p.x[2] != cp.x[2]))
				return p->p;
			else if (j == i)
				return p->next->p;
			else
				p = p->next;
		}
	}

	inline Particle& operator()(size_t i)
	{
		assert(i < length);
		ParticleNode * p = header;
		int j;
		for ( j = 0; j < length; j++)
		{
			if (j == i)
				return p->p;
			else
				p = p->next;
		}		
	}

	inline Particle& operator()(size_t i) const
	{
		assert(i < length);
		ParticleNode * p = header;
		int j;
		for (j = 0; j < length; j++)
		{
			if (j == i)
				return p->p;
			else
				p = p->next;
		}
	}

	void addParticle(Particle p)
	{
		ParticleNode* pn;
		pn->p = p;
		pn->next = header;
		header = pn;
		length++;
	}

	bool removeParticle(ParticleNode *p)
	{
		ParticleNode* pn = header;
		ParticleNode* pnp = header;
		int i;
		for ( i = 0; i < length; i++)
		{
			if (pn == p)
			{
				pnp->next = pn->next;
				delete p;
				length--;
				return true;
			}
			else
			{
				pnp = pn;
				pn = pn->next;
			}
		}
		if (i == length)
			return false;
	}

	bool removeParticle(size_t p)
	{
		ParticleNode* pn = header;
		ParticleNode* pnp = header;
		int i;
		for (i = 0; i < length; i++)
		{
			if (i == p)
			{
				pnp->next = pn->next;
				delete pn;
				length--;
				return true;
			}
			else
			{
				pnp = pn;
				pn = pn->next;
			}
		}
		if (i == length)
			return false;
	}

	bool insertAfter(ParticleNode *sp, ParticleNode *tp)
	{
		ParticleNode* pn = header;
		int i;
		for (i = 0; i < length; i++)
		{
			if (pn == sp)
			{
				tp->next = sp->next;
				sp->next = tp;
				length++;
				return true;
			}
			else
			{
				pn = pn->next;
			}
		}

		if (i == length)
			return false;
	}

};