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

#ifndef PE_APP_H
#define PE_APP_H

#ifdef __cplusplus
extern "C" {
#endif

/************************************************************************
 ** Pragmas
 *************************************************************************/

/************************************************************************
 ** Includes
 *************************************************************************/

#include "cfe.h"

#include "pe_platform_cfg.h"
#include "pe_mission_cfg.h"
#include "pe_perfids.h"
#include "pe_msgids.h"
#include "pe_msg.h"
#include "pe_events.h"
#include "pe_tbldefs.h"
#include "px4_msgs.h"
#include "px4lib.h"

#include "math/Vector1F.hpp"
#include "math/Vector3F.hpp"
#include "math/Vector10F.hpp"
#include "math/Matrix10F10.hpp"
#include "math/Matrix10F3.hpp"
#include "math/Matrix3F3.hpp"
#include "math/LowPass.hpp"
#include "math/LowPassVector10F.hpp"
#include "math/Stats1F.hpp"
#include "math/Stats6F.hpp"
#include "math/Quaternion.hpp"
#include "math/Euler.hpp"
#include "math/Dcm.hpp"
#include "geo/geo.h"

#include <poll.h>
#include <math.h>

/************************************************************************
 ** Local Defines
 *************************************************************************/

#define PE_BETA_TABLE_SIZE	(7)

/************************************************************************
 ** Local Structure Definitions
 *************************************************************************/

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
	float  LAND_VXY_STDDEV;// TODO TODO TODO XXX Check what this is supposed to be
	float  PN_P_NOISE_DENSITY;
	float  PN_V_NOISE_DENSITY;
	float  PN_B_NOISE_DENSITY;
	float  PN_T_NOISE_DENSITY;
	float  T_MAX_GRADE;
	int32  FAKE_ORIGIN;
	float  INIT_ORIGIN_LAT;
	float  INIT_ORIGIN_LON;
} PE_Params_t;

enum {
	X_x  = 0,
	X_y  = 1,
	X_z  = 2,
	X_vx = 3,
	X_vy = 4,
	X_vz = 5,
	X_bx = 6,
	X_by = 7,
	X_bz = 8,
	X_tz = 9,
	n_x  = 10
};

enum {
	U_ax = 0,
	U_ay = 1,
	U_az = 2,
	n_u  = 3
};

enum {
	Y_baro_z = 0,
	n_y_baro = 1
};

enum {
	Y_gps_x  = 0,
	Y_gps_y  = 1,
	Y_gps_z  = 2,
	Y_gps_vx = 3,
	Y_gps_vy = 4,
	Y_gps_vz = 5,
	n_y_gps  = 6
};

enum {
	Y_land_vx  = 0,
	Y_land_vy  = 1,
	Y_land_agl = 2,
	n_y_land   = 3
};

/* Enums for other sensors would go here */


/**
 **  \brief PE Application Class
 */
class PE
{
public:
    PE();
    ~PE();

    /** \brief Constants */
    float 	DELAY_MAX;
    float 	HIST_STEP;
    float 	BIAS_MAX;
    size_t 	HIST_LEN = 10;
    size_t 	N_DIST_SUBS;
    float  	BETA_TABLE[PE_BETA_TABLE_SIZE];
    uint32 	EST_STDDEV_XY_VALID; // 2.0 m
    uint32 	EST_STDDEV_Z_VALID; // 2.0 m
    uint32 	EST_STDDEV_TZ_VALID; // 2.0 m
    float 	P_MAX; // max allowed value in state covariance
    float 	LAND_RATE; // rate of land detector correction

    /**\brief Scheduling Pipe ID */
    CFE_SB_PipeId_t SchPipeId;

    /** \brief Command Pipe ID */
    CFE_SB_PipeId_t CmdPipeId;

    /* Task-related */

    /** \brief Task Run Status */
    uint32 uiRunStatus;

    /* Config table-related */

    /** \brief Config Table Handle */
    CFE_TBL_Handle_t ConfigTblHdl;

    /** \brief Config Table Pointer */
    PE_ConfigTbl_t* ConfigTblPtr;

    int32 ConfigMutex;

    /** \brief Ingest Data */
    PX4_VehicleGpsPositionMsg_t m_VehicleGpsPositionMsg;
    PX4_VehicleStatusMsg_t m_VehicleStatusMsg; //todo verify needed
    PX4_VehicleLandDetectedMsg_t m_VehicleLandDetectedMsg;
    PX4_ActuatorArmedMsg_t m_ActuatorArmedMsg;
    PX4_VehicleAttitudeMsg_t m_VehicleAttitudeMsg;
    PX4_VehicleControlModeMsg_t m_VehicleControlModeMsg; //todo verify needed
    PX4_SensorCombinedMsg_t m_SensorCombinedMsg;
    PX4_VehicleAttitudeSetpointMsg_t m_VehicleAttitudeSetpointMsg;
    PX4_ManualControlSetpointMsg_t m_ManualControlSetpointMsg; //todo verify needed

    /** \brief Output Data published at the end of cycle */
    PX4_VehicleLocalPositionMsg_t m_VehicleLocalPositionMsg;
    PX4_EstimatorStatusMsg_t m_EstimatorStatusMsg;
    PX4_VehicleGlobalPositionMsg_t m_VehicleGlobalPositionMsg;
    PX4_Ekf2InnovationsMsg_t m_Ekf2InnovationsMsg;

    bool m_ParamsUpdated;

    // Sensor stats
    Stats1F m_BaroStats;
    Stats6F m_GpsStats;
    uint16 m_LandCount;

    bool m_XyEstValid;
    bool m_ZEstValid;
    bool m_TzEstValid;

    // map
    struct map_projection_reference_s m_MapRef;

	// low pass
	LowPassVector10F m_XLowPass;
	LowPass m_AglLowPass;

	// delay blocks
	//BlockDelay<float, n_x, 1, HIST_LEN> m_XDelay;
	//BlockDelay<uint64_t, 1, 1, HIST_LEN> m_TDelay;

	// Time
    pollfd m_Polls[3]; // don't think needed
	uint64 m_Timestamp;
	uint64 m_TimestampLastBaro;
	uint64 m_TimeLastBaro;
	uint64 m_TimeLastGps;
	uint64 m_TimeLastLand;
	bool   m_BaroTimeout;
	bool   m_GpsTimeout;
	bool   m_LandTimeout;

	// reference altitudes
	float m_AltOrigin;
	bool  m_AltOriginInitialized;
	float m_BaroAltOrigin;
	float m_GpsAltOrigin;

	// status
	bool m_ReceivedGps;
	bool m_LastArmedState;

	// masks
	uint8 m_SensorTimeout;
	uint8 m_SensorFault;
	bool m_EstimatorInitialized;

	// state space
	math::Vector10F  m_StateVec; // state vector
	math::Vector3F  m_InputVec; // input vector
	math::Matrix10F10  m_StateCov; // state covariance matrix

	math::Matrix3F3 _R_att;
	math::Vector3F m_Euler;

	math::Matrix10F10 m_DynamicsMat; // dynamics matrix
	math::Matrix10F3 m_InputMat; // input matrix
	math::Matrix3F3 m_InputCov; // input covariance
	math::Matrix10F10 m_NoiseCov; // process noise covariance

	PE_Params_t m_Params;

    /** \brief Housekeeping Telemetry for downlink */
    PE_HkTlm_t HkTlm;

    /************************************************************************/
    /** \brief Position Estimator (PE) application entry point
     **
     **  \par Description
     **       Position Estimator Task application entry point.  This function
     **       performs app initialization, then waits for the cFE ES Startup
     **       Sync, then executes the main processing loop.
     **
     **  \par Assumptions, External Events, and Notes:
     **       If there is an unrecoverable failure during initialization the
     **       main loop is never executed and the application will exit.
     **
     *************************************************************************/
    void AppMain(void);

    /************************************************************************/
    /** \brief Initialize the Position Estimator (PE) application
     **
     **  \par Description
     **       Position Estimator application initialization routine. This
     **       function performs all the required startup steps to
     **       initialize (or restore from CDS) PE data structures and get
     **       the application registered with the cFE services so it can
     **       begin to receive command messages and send events.
     **
     **  \par Assumptions, External Events, and Notes:
     **       None
     **
     **  \returns
     **  \retcode #CFE_SUCCESS  \retdesc \copydoc CFE_SUCCESS    \endcode
     **  \retstmt Return codes from #CFE_ES_RegisterApp          \endcode
     **  \retstmt Return codes from #PE_InitEvent               \endcode
     **  \retstmt Return codes from #PE_InitPipe                \endcode
     **  \retstmt Return codes from #PE_InitData                \endcode
     **  \retstmt Return codes from #PE_InitConfigTbl           \endcode
     **  \retstmt Return codes from #OS_TaskInstallDeleteHandler \endcode
     **  \endreturns
     **
     *************************************************************************/
    int32 InitApp(void);

    /************************************************************************/
    /** \brief Initialize Event Services and Event tables
     **
     **  \par Description
     **       This function performs the steps required to setup
     **       cFE Event Services for use by the PE application.
     **
     **  \par Assumptions, External Events, and Notes:
     **       None
     **
     **  \returns
     **  \retcode #CFE_SUCCESS  \retdesc \copydoc CFE_SUCCESS \endcode
     **  \retstmt Return codes from #CFE_EVS_Register  \endcode
     **  \endreturns
     **
     *************************************************************************/
    int32 InitEvent(void);

    /************************************************************************/
    /** \brief Initialize global variables used by PE application
     **
     **  \par Description
     **       This function performs the steps required to initialize
     **       the PE application data.
     **
     **  \par Assumptions, External Events, and Notes:
     **       None
     **
     *************************************************************************/
    void InitData(void);

    /************************************************************************/
    /** \brief Initialize message pipes
     **
     **  \par Description
     **       This function performs the steps required to setup
     **       initialize the cFE Software Bus message pipes and subscribe to
     **       messages for the PE application.
     **
     **  \par Assumptions, External Events, and Notes:
     **       None
     **
     **  \returns
     **  \retcode #CFE_SUCCESS  \retdesc \copydoc CFE_SUCCESS \endcode
     **  \retstmt Return codes from #CFE_SB_CreatePipe        \endcode
     **  \retstmt Return codes from #CFE_SB_SubscribeEx       \endcode
     **  \retstmt Return codes from #CFE_SB_Subscribe         \endcode
     **  \endreturns
     **
     *************************************************************************/
    int32 InitPipe(void);

    /************************************************************************/
    /** \brief Receive and process messages from the scheduler pipe.
     **
     **  \par Description
     **       This function receives and processes messages
     **       for the PE application from the SCH pipe.  This function
     **       will pend for the type defined by iBlocking, allowing
     **       it to wait for messages, i.e. wakeup messages from scheduler.
     **
     **  \par Assumptions, External Events, and Notes:
     **       None
     **
     **  \param [in]   iBlocking    A #CFE_SB_PEND_FOREVER, #CFE_SB_POLL or
     **                             millisecond timeout
     **
     **  \returns
     **  \retcode #CFE_SUCCESS  \retdesc \copydoc CFE_SUCCESS \endcode
     **  \retstmt Return codes from #CFE_SB_RcvMsg            \endcode
     **  \endreturns
     **
     *************************************************************************/
    int32 RcvSchPipeMsg(int32 iBlocking);

    /************************************************************************/
    /** \brief Position Estimator Task incoming command processing
     **
     **  \par Description
     **       This function processes incoming commands subscribed
     **       by PE application
     **
     **  \par Assumptions, External Events, and Notes:
     **       None
     **
     *************************************************************************/
    void ProcessCmdPipe(void);

    /************************************************************************/
    /** \brief Position Estimator Task application commands
     **
     **  \par Description
     **       This function processes command messages
     **       specific to the PE application
     **
     **  \par Assumptions, External Events, and Notes:
     **       None
     **
     **  \param [in]   MsgPtr       A #CFE_SB_Msg_t pointer that
     **                             references the software bus message
     **
     *************************************************************************/
    void ProcessAppCmds(CFE_SB_Msg_t* MsgPtr);

    /************************************************************************/
    /** \brief Sends PE housekeeping message
     **
     **  \par Description
     **       This function sends the housekeeping message
     **
     **  \par Assumptions, External Events, and Notes:
     **       None
     **
     *************************************************************************/
    void ReportHousekeeping(void);
    /************************************************************************/
    /** \brief Sends the VehicleLocalPositionMsg message.
     **
     **  \par Description
     **       This function publishes the VehicleLocalPositionMsg message containing
     **       <TODO>
     **
     **  \par Assumptions, External Events, and Notes:
     **       None
     **
     *************************************************************************/
    void SendVehicleLocalPositionMsg(void);
    /************************************************************************/
    /** \brief Sends the EstimatorStatusMsg message.
     **
     **  \par Description
     **       This function publishes the EstimatorStatusMsg message containing
     **       <TODO>
     **
     **  \par Assumptions, External Events, and Notes:
     **       None
     **
     *************************************************************************/
    void SendEstimatorStatusMsg(void);
    /************************************************************************/
    /** \brief Sends the VehicleGlobalPositionMsg message.
     **
     **  \par Description
     **       This function publishes the VehicleGlobalPositionMsg message containing
     **       <TODO>
     **
     **  \par Assumptions, External Events, and Notes:
     **       None
     **
     *************************************************************************/
    void SendVehicleGlobalPositionMsg(void);
    /************************************************************************/
    /** \brief Sends the Ekf2InnovationsMsg message.
     **
     **  \par Description
     **       This function publishes the Ekf2InnovationsMsg message containing
     **       <TODO>
     **
     **  \par Assumptions, External Events, and Notes:
     **       None
     **
     *************************************************************************/
    void SendEkf2InnovationsMsg(void);
    /************************************************************************/
    /** \brief Verify Command Length
     **
     **  \par Description
     **       This function verifies the command message length.
     **
     **  \par Assumptions, External Events, and Notes:
     **       None
     **
     **  \param [in]   MsgPtr        A #CFE_SB_Msg_t pointer that
     **                              references the software bus message
     **  \param [in]   usExpectedLen The expected length of the message
     **
     **  \returns
     **  TRUE if the message length matches expectations, FALSE if it does not.
     **  \endreturns
     **
     *************************************************************************/
    boolean VerifyCmdLength(CFE_SB_Msg_t* MsgPtr, uint16 usExpectedLen);

    /************************************************************************/
    /** \brief Validate PE configuration table
    **
    **  \par Description
    **       This function validates PE's configuration table
    **
    **  \par Assumptions, External Events, and Notes:
    **       None
    **
    **  \param [in]   ConfigTblPtr    A pointer to the table to validate.
    **
    **  \returns
    **  \retcode #CFE_SUCCESS  \retdesc \copydoc CFE_SUCCESS  \endcode
    **  \endreturns
    **
    *************************************************************************/
    static int32  ValidateConfigTbl(void*);

private:
    /************************************************************************/
    /** \brief Initialize the PE configuration tables.
    **
    **  \par Description
    **       This function initializes PE's configuration tables.  This
    **       includes <TODO>.
    **
    **  \par Assumptions, External Events, and Notes:
    **       None
    **
    **  \returns
    **  \retcode #CFE_SUCCESS  \retdesc \copydoc CFE_SUCCESS  \endcode
    **  \retstmt Return codes from #CFE_TBL_Register          \endcode
    **  \retstmt Return codes from #CFE_TBL_Load              \endcode
    **  \retstmt Return codes from #PE_AcquireConfigPointers \endcode
    **  \endreturns
    **
    *************************************************************************/
    int32  InitConfigTbl(void);

    /************************************************************************/
    /** \brief Obtain PE configuration tables data pointers.
    **
    **  \par Description
    **       This function manages the configuration tables
    **       and obtains a pointer to their data.
    **
    **  \par Assumptions, External Events, and Notes:
    **       None
    **
    **  \returns
    **  \retcode #CFE_SUCCESS  \retdesc \copydoc CFE_SUCCESS  \endcode
    **  \endreturns
    **
    *************************************************************************/
    int32  AcquireConfigPointers(void);



private:
    //Vector<n_x> dynamics(const Vector<float, n_x> &x,const Vector<float, n_u> &u);
	void initStateCov();
	void InitStateSpace();
	void updateStateSpace();
	void updateStateSpaceParams();

	// predict the next state
	void predict();

	// baro
	int32  baroMeasure(math::Vector1F &y);
	void baroCorrect();
	void baroInit();
	void baroCheckTimeout();

	// gps
	int  gpsMeasure(math::Vector6F &y);
	void gpsCorrect();
	void gpsInit();
	void gpsCheckTimeout();

	// land
	int  landMeasure(math::Vector3F &y);
	void landCorrect();
	void landInit();
	void landCheckTimeout();

	// timeouts
	void CheckTimeouts();

	// misc
	//inline float agl() { return _x(X_tz) - _x(X_z); }
	bool landed();
	int getDelayPeriods(float delay, uint8 *periods);

	void UpdateLocalParams();

	void Update();
	void Predict(float dt);

	math::Vector10F dynamics(const math::Vector10F &x, const math::Vector3F &u);






};

#ifdef __cplusplus
}
#endif 

#endif /* PE_APP_H */

/************************/
/*  End of File Comment */
/************************/
