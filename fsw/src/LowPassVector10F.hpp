#ifndef LOWPASSVECTOR10F_H
#define LOWPASSVECTOR10F_H

#include "cfe.h"
#include "Vector10F.hpp"
#include "Matrix10F10.hpp"

#define SIZE    (10)

class LowPassVector10F
{
public:
	LowPassVector10F();
	~LowPassVector10F();
	math::Vector10F Update(math::Matrix10F10 input, float dt, float cutoff);

private:
	bool m_Initialized;
	math::Vector10F m_State;
};

#endif
