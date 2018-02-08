#ifndef LOWPASS_H
#define LOWPASS_H

#include "cfe.h"

class LowPass
{
public:
	LowPass();
	~LowPass();
	float Update(float value, float dt, float cutoff);

private:
	bool m_Initialized;
	float m_State;
};

#endif
