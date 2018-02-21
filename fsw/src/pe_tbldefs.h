/****************************************************************************
*
*   Copyright (c) 2017 Windhover Labs, L.L.C. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in
*    the documentation and/or other materials provided with the
*    distribution.
* 3. Neither the name Windhover Labs nor the names of its
*    contributors may be used to endorse or promote products derived
*    from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
* COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/

#ifndef PE_TBLDEFS_H
#define PE_TBLDEFS_H

/************************************************************************
** Pragmas
*************************************************************************/

/************************************************************************
** Includes
*************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif

#include "cfe.h"
#include "pe_platform_cfg.h"

/************************************************************************
** Local Defines
*************************************************************************/

/**
 * \brief Defines the table identification name used for the
 * configuration table registration.
 */
#define PE_CONFIG_TABLENAME ("CONFIG_TBL")


/** \brief Params table */
typedef struct
{
	int32  FUSION;

	float  VXY_PUB_THRESH;

	float  Z_PUB_THRESH;

	float  ACCEL_XY_STDDEV;

	float  ACCEL_Z_STDDEV;

	float  BARO_STDDEV;

	float  GPS_DELAY;

	float  GPS_XY_STDDEV;

	float  GPS_Z_STDDEV;

	float  GPS_VXY_STDDEV;

	float  GPS_VZ_STDDEV;

	float  GPS_EPH_MAX;

	float  GPS_EPV_MAX;

	float  LAND_Z_STDDEV;

	float  LAND_VXY_STDDEV;

	float  PN_P_NOISE_DENSITY;

	float  PN_V_NOISE_DENSITY;

	float  PN_B_NOISE_DENSITY;

	float  PN_T_NOISE_DENSITY;

	float  T_MAX_GRADE;

	int32  FAKE_ORIGIN;

	float  INIT_ORIGIN_LAT;

	float  INIT_ORIGIN_LON;

} PE_ConfigTbl_t;


#ifdef __cplusplus
}
#endif

#endif /* PE_TBLDEFS_H */

/************************/
/*  End of File Comment */
/************************/
