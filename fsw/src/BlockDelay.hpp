#ifndef BLOCK_DELAY_HPP
#define BLOCK_DELAY_HPP

#include "cfe.h"
#include <math/Matrix10F1.hpp>

namespace delay
{

class BlockDelay10F1LEN10
{
private:
    static const int LEN = 10;
    uint8 m_Index;
    int m_Delay;
    math::Matrix10F1 m_Storage[LEN];
    

public:
    BlockDelay10F1LEN10();
    ~BlockDelay10F1LEN10();
    math::Matrix10F1 Update(const math::Matrix10F1 &input);
    void Initialize(void);
    math::Matrix10F1 Get(uint8 delay);

protected:

};


class BlockDelayUINT64LEN10
{
private:
    static const int LEN = 10;
    uint8 m_Index;
    int m_Delay;
    uint64 m_Storage[LEN];
    

public:
    BlockDelayUINT64LEN10();
    ~BlockDelayUINT64LEN10();
    uint64 Update(const uint64 &input);
    void Initialize(void);
    uint64 Get(uint8 delay);

protected:

};

}

#endif
