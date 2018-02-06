/************************************************************************
** Includes
*************************************************************************/
#include <string.h>

#include "cfe.h"

#include "pe_app.h"
#include "pe_msg.h"
#include "pe_version.h"
#include "lib/px4lib.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* Instantiate the application object.                             */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
PE oPE;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* Default constructor.                                            */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
PE::PE()
{

}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* Destructor constructor.                                         */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
PE::~PE()
{

}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* Initialize event tables.                                        */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int32 PE::InitEvent()
{
    int32  iStatus=CFE_SUCCESS;

    /* Register the table with CFE */
    iStatus = CFE_EVS_Register(0, 0, CFE_EVS_BINARY_FILTER);
    if (iStatus != CFE_SUCCESS)
    {
        (void) CFE_ES_WriteToSysLog("PE - Failed to register with EVS (0x%08lX)\n", iStatus);
    }

    return iStatus;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* Initialize Message Pipes                                        */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int32 PE::InitPipe()
{
    int32  iStatus=CFE_SUCCESS;

    /* Init schedule pipe and subscribe to wakeup messages */
    iStatus = CFE_SB_CreatePipe(&SchPipeId,
    		PE_SCH_PIPE_DEPTH,
			PE_SCH_PIPE_NAME);
    if (iStatus == CFE_SUCCESS)
    {
        iStatus = CFE_SB_SubscribeEx(PE_WAKEUP_MID, SchPipeId, CFE_SB_Default_Qos, PE_WAKEUP_MID_MAX_MSG_COUNT);
        if (iStatus != CFE_SUCCESS)
        {
            (void) CFE_EVS_SendEvent(PE_SUBSCRIBE_ERR_EID, CFE_EVS_ERROR,
            		"Sch Pipe failed to subscribe to PE_WAKEUP_MID. (0x%08lX)",
                    iStatus);
            goto PE_InitPipe_Exit_Tag;
        }

        iStatus = CFE_SB_SubscribeEx(PE_SEND_HK_MID, SchPipeId, CFE_SB_Default_Qos, PE_SEND_HK_MID_MAX_MSG_COUNT);
        if (iStatus != CFE_SUCCESS)
        {
            (void) CFE_EVS_SendEvent(PE_SUBSCRIBE_ERR_EID, CFE_EVS_ERROR,
					 "CMD Pipe failed to subscribe to PE_SEND_HK_MID. (0x%08X)",
					 (unsigned int)iStatus);
            goto PE_InitPipe_Exit_Tag;
        }
        iStatus = CFE_SB_SubscribeEx(PX4_VEHICLE_GPS_POSITION_MID, SchPipeId, CFE_SB_Default_Qos, 1);
        if (iStatus != CFE_SUCCESS)
        {
            (void) CFE_EVS_SendEvent(PE_SUBSCRIBE_ERR_EID, CFE_EVS_ERROR,
					 "CMD Pipe failed to subscribe to PX4_VEHICLE_GPS_POSITION_MID. (0x%08lX)",
					 iStatus);
            goto PE_InitPipe_Exit_Tag;
        }
        iStatus = CFE_SB_SubscribeEx(PX4_VEHICLE_STATUS_MID, SchPipeId, CFE_SB_Default_Qos, 1);
        if (iStatus != CFE_SUCCESS)
        {
            (void) CFE_EVS_SendEvent(PE_SUBSCRIBE_ERR_EID, CFE_EVS_ERROR,
					 "CMD Pipe failed to subscribe to PX4_VEHICLE_STATUS_MID. (0x%08lX)",
					 iStatus);
            goto PE_InitPipe_Exit_Tag;
        }
        iStatus = CFE_SB_SubscribeEx(PX4_VEHICLE_LAND_DETECTED_MID, SchPipeId, CFE_SB_Default_Qos, 1);
        if (iStatus != CFE_SUCCESS)
        {
            (void) CFE_EVS_SendEvent(PE_SUBSCRIBE_ERR_EID, CFE_EVS_ERROR,
					 "CMD Pipe failed to subscribe to PX4_VEHICLE_LAND_DETECTED_MID. (0x%08lX)",
					 iStatus);
            goto PE_InitPipe_Exit_Tag;
        }
        iStatus = CFE_SB_SubscribeEx(PX4_ACTUATOR_ARMED_MID, SchPipeId, CFE_SB_Default_Qos, 1);
        if (iStatus != CFE_SUCCESS)
        {
            (void) CFE_EVS_SendEvent(PE_SUBSCRIBE_ERR_EID, CFE_EVS_ERROR,
					 "CMD Pipe failed to subscribe to PX4_ACTUATOR_ARMED_MID. (0x%08lX)",
					 iStatus);
            goto PE_InitPipe_Exit_Tag;
        }
        iStatus = CFE_SB_SubscribeEx(PX4_VEHICLE_ATTITUDE_MID, SchPipeId, CFE_SB_Default_Qos, 1);
        if (iStatus != CFE_SUCCESS)
        {
            (void) CFE_EVS_SendEvent(PE_SUBSCRIBE_ERR_EID, CFE_EVS_ERROR,
					 "CMD Pipe failed to subscribe to PX4_VEHICLE_ATTITUDE_MID. (0x%08lX)",
					 iStatus);
            goto PE_InitPipe_Exit_Tag;
        }
        iStatus = CFE_SB_SubscribeEx(PX4_VEHICLE_CONTROL_MODE_MID, SchPipeId, CFE_SB_Default_Qos, 1);
        if (iStatus != CFE_SUCCESS)
        {
            (void) CFE_EVS_SendEvent(PE_SUBSCRIBE_ERR_EID, CFE_EVS_ERROR,
					 "CMD Pipe failed to subscribe to PX4_VEHICLE_CONTROL_MODE_MID. (0x%08lX)",
					 iStatus);
            goto PE_InitPipe_Exit_Tag;
        }
        iStatus = CFE_SB_SubscribeEx(PX4_SENSOR_COMBINED_MID, SchPipeId, CFE_SB_Default_Qos, 1);
        if (iStatus != CFE_SUCCESS)
        {
            (void) CFE_EVS_SendEvent(PE_SUBSCRIBE_ERR_EID, CFE_EVS_ERROR,
					 "CMD Pipe failed to subscribe to PX4_SENSOR_COMBINED_MID. (0x%08lX)",
					 iStatus);
            goto PE_InitPipe_Exit_Tag;
        }
        iStatus = CFE_SB_SubscribeEx(PX4_VEHICLE_ATTITUDE_SETPOINT_MID, SchPipeId, CFE_SB_Default_Qos, 1);
        if (iStatus != CFE_SUCCESS)
        {
            (void) CFE_EVS_SendEvent(PE_SUBSCRIBE_ERR_EID, CFE_EVS_ERROR,
					 "CMD Pipe failed to subscribe to PX4_VEHICLE_ATTITUDE_SETPOINT_MID. (0x%08lX)",
					 iStatus);
            goto PE_InitPipe_Exit_Tag;
        }
        iStatus = CFE_SB_SubscribeEx(PX4_MANUAL_CONTROL_SETPOINT_MID, SchPipeId, CFE_SB_Default_Qos, 1);
        if (iStatus != CFE_SUCCESS)
        {
            (void) CFE_EVS_SendEvent(PE_SUBSCRIBE_ERR_EID, CFE_EVS_ERROR,
					 "CMD Pipe failed to subscribe to PX4_MANUAL_CONTROL_SETPOINT_MID. (0x%08lX)",
					 iStatus);
            goto PE_InitPipe_Exit_Tag;
        }
        iStatus = CFE_SB_SubscribeEx(PX4_DISTANCE_SENSOR_MID, SchPipeId, CFE_SB_Default_Qos, 1);
        if (iStatus != CFE_SUCCESS)
        {
            (void) CFE_EVS_SendEvent(PE_SUBSCRIBE_ERR_EID, CFE_EVS_ERROR,
					 "CMD Pipe failed to subscribe to PX4_DISTANCE_SENSOR_MID. (0x%08lX)",
					 iStatus);
            goto PE_InitPipe_Exit_Tag;
        }
    }
    else
    {
        (void) CFE_EVS_SendEvent(PE_PIPE_INIT_ERR_EID, CFE_EVS_ERROR,
			 "Failed to create SCH pipe (0x%08lX)",
			 iStatus);
        goto PE_InitPipe_Exit_Tag;
    }

    /* Init command pipe and subscribe to command messages */
    iStatus = CFE_SB_CreatePipe(&CmdPipeId,
    		PE_CMD_PIPE_DEPTH,
			PE_CMD_PIPE_NAME);
    if (iStatus == CFE_SUCCESS)
    {
        /* Subscribe to command messages */
        iStatus = CFE_SB_Subscribe(PE_CMD_MID, CmdPipeId);

        if (iStatus != CFE_SUCCESS)
        {
            (void) CFE_EVS_SendEvent(PE_SUBSCRIBE_ERR_EID, CFE_EVS_ERROR,
				 "CMD Pipe failed to subscribe to PE_CMD_MID. (0x%08lX)",
				 iStatus);
            goto PE_InitPipe_Exit_Tag;
        }
    }
    else
    {
        (void) CFE_EVS_SendEvent(PE_PIPE_INIT_ERR_EID, CFE_EVS_ERROR,
			 "Failed to create CMD pipe (0x%08lX)",
			 iStatus);
        goto PE_InitPipe_Exit_Tag;
    }

PE_InitPipe_Exit_Tag:
    return iStatus;
}
    

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* Initialize Global Variables                                     */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void PE::InitData()
{
    /* Init housekeeping message. */
    CFE_SB_InitMsg(&HkTlm,
    		PE_HK_TLM_MID, sizeof(HkTlm), TRUE);

    /* Init output messages */
	CFE_SB_InitMsg(&mVehicleLocalPositionMsg,
			PX4_VEHICLE_LOCAL_POSITION_MID, sizeof(PX4_VehicleLocalPositionMsg_t), TRUE);
	CFE_SB_InitMsg(&mEstimatorStatusMsg,
			PX4_ESTIMATOR_STATUS_MID, sizeof(PX4_EstimatorStatusMsg_t), TRUE);
	CFE_SB_InitMsg(&mVehicleGlobalPositionMsg,
			PX4_VEHICLE_GLOBAL_POSITION_MID, sizeof(PX4_VehicleGlobalPositionMsg_t), TRUE);
	CFE_SB_InitMsg(&mEkf2InnovationsMsg,
			PX4_EKF2_INNOVATIONS_MID, sizeof(PX4_Ekf2InnovationsMsg_t), TRUE);

	/* Set constants */
	DELAY_MAX = 0.5f; // seconds
	HIST_STEP = 0.05f; // 20 hz
	BIAS_MAX = 1e-1f;
	HIST_LEN = 10; // DELAY_MAX / HIST_STEP;
	N_DIST_SUBS = 4;
	BETA_TABLE[0] = 0;
	BETA_TABLE[1] = 8.82050518214;
	BETA_TABLE[2] = 12.094592431;
	BETA_TABLE[3] = 13.9876612368;
	BETA_TABLE[4] = 16.0875642296;
	BETA_TABLE[5] = 17.8797700658;
	BETA_TABLE[6] = 19.6465647819;
	EST_STDDEV_XY_VALID = 2.0;
	EST_STDDEV_Z_VALID = 2.0;
	EST_STDDEV_TZ_VALID = 2.0;
	P_MAX = 1.0e6f;
	LAND_RATE = 10.0f;

	/* Map projection */
	// todo

	/* Block stats? todo */
	//mBaroStats(this, "");
	//mGpsStats(this, "");

	/* Block delays? todo */
	//mXDelay(this, "");
	//mTDelay(this, "");

	/*  */
	memset(&mPolls, 0, sizeof(mPolls));
	mTimeStamp = PX4LIB_GetPX4TimeUs();
	mTimeStampLastBaro = PX4LIB_GetPX4TimeUs();
	mTimeLastBaro = 0;
	mTimeLastGps = 0;
	mTimeLastLand = 0;

	// reference altitudes
	mAltOrigin = 0.0;
	mAltOriginInitialized = FALSE;
	mBaroAltOrigin = 0.0;
	mGpsAltOrigin = 0.0;

	// status
	mReceivedGps = FALSE;
	mLastArmedState = FALSE;

	// masks todo
	mSensorTimeout = 255;
	mSensorFault = 0;
	mEstimatorInitialized = 0;

	// state space
	//mStateVec.Zero();
	//mInputVec.Zero();
	//Matrix<float, n_x, n_x>  mStateCov; // state covariance matrix

	//matrix::Dcm<float> _R_att;
	//Vector3f mEuler;

	//Matrix<float, n_x, n_x>  mDynamicsMat; // dynamics matrix
	//Matrix<float, n_x, n_u>  mInputMat; // input matrix
	//Matrix<float, n_u, n_u>  mInputCov; // input covariance
	//Matrix<float, n_x, n_x>  mNoiseCov; // process noise covariance

	/* Map */
	//mMap_ref.init_done = false;

	/*  */


}

void PE::initStateCov()
{
//	mStateCov.setZero();
//	// initialize to twice valid condition
//	mStateCov(X_x, X_x) = 2 * EST_STDDEV_XY_VALID * EST_STDDEV_XY_VALID;
//	mStateCov(X_y, X_y) = 2 * EST_STDDEV_XY_VALID * EST_STDDEV_XY_VALID;
//	mStateCov(X_z, X_z) = 2 * EST_STDDEV_Z_VALID * EST_STDDEV_Z_VALID;
//	mStateCov(X_vx, X_vx) = 2 * tbl->vxy_pub_thresh * tbl->vxy_pub_thresh;
//	mStateCov(X_vy, X_vy) = 2 * tbl->vxy_pub_thresh * tbl->vxy_pub_thresh;
//	// use vxy thresh for vz init as well
//	mStateCov(X_vz, X_vz) = 2 * tbl->vxy_pub_thresh * tbl->vxy_pub_thresh;
//	// initialize bias uncertainty to small values to keep them stable
//	mStateCov(X_bx, X_bx) = 1e-6;
//	mStateCov(X_by, X_by) = 1e-6;
//	mStateCov(X_bz, X_bz) = 1e-6;
//	mStateCov(X_tz, X_tz) = 2 * EST_STDDEV_TZ_VALID * EST_STDDEV_TZ_VALID;
}

void PE::InitStateSpace()
{
	initStateCov();

//	// dynamics matrix
//	mDynamicsMat.Zero();
//	// derivative of position is velocity
//	mDynamicsMat(X_x, X_vx) = 1;
//	mDynamicsMat(X_y, X_vy) = 1;
//	mDynamicsMat(X_z, X_vz) = 1;
//
//	// input matrix
//	mInputMat.Zero();
//	mInputMat(X_vx, U_ax) = 1;
//	mInputMat(X_vy, U_ay) = 1;
//	mInputMat(X_vz, U_az) = 1;

	// update components that depend on current state
	updateStateSpace();
	updateStateSpaceParams();
}

void PE::updateStateSpace()
{
	// derivative of velocity is accelerometer acceleration
	// (in input matrix) - bias (in body frame)
//	mDynamicsMat(X_vx, X_bx) = -_R_att(0, 0);
//	mDynamicsMat(X_vx, X_by) = -_R_att(0, 1);
//	mDynamicsMat(X_vx, X_bz) = -_R_att(0, 2);
//
//	mDynamicsMat(X_vy, X_bx) = -_R_att(1, 0);
//	mDynamicsMat(X_vy, X_by) = -_R_att(1, 1);
//	mDynamicsMat(X_vy, X_bz) = -_R_att(1, 2);
//
//	mDynamicsMat(X_vz, X_bx) = -_R_att(2, 0);
//	mDynamicsMat(X_vz, X_by) = -_R_att(2, 1);
//	mDynamicsMat(X_vz, X_bz) = -_R_att(2, 2);
}

void PE::updateStateSpaceParams()
{
//	// input noise covariance matrix
//	mInputCov.setZero();
//	mInputCov(U_ax, U_ax) = _accel_xy_stddev.get() * _accel_xy_stddev.get();
//	mInputCov(U_ay, U_ay) = _accel_xy_stddev.get() * _accel_xy_stddev.get();
//	mInputCov(U_az, U_az) = _accel_z_stddev.get() * _accel_z_stddev.get();
//
//	// process noise power matrix
//	mNoiseCov.setZero();
//	float pn_p_sq = _pn_p_noise_density.get() * _pn_p_noise_density.get();
//	float pn_v_sq = _pn_v_noise_density.get() * _pn_v_noise_density.get();
//	mNoiseCov(X_x, X_x) = pn_p_sq;
//	mNoiseCov(X_y, X_y) = pn_p_sq;
//	mNoiseCov(X_z, X_z) = pn_p_sq;
//	mNoiseCov(X_vx, X_vx) = pn_v_sq;
//	mNoiseCov(X_vy, X_vy) = pn_v_sq;
//	mNoiseCov(X_vz, X_vz) = pn_v_sq;
//
//	// technically, the noise is in the body frame,
//	// but the components are all the same, so
//	// ignoring for now
//	float pn_b_sq = _pn_b_noise_density.get() * _pn_b_noise_density.get();
//	mNoiseCov(X_bx, X_bx) = pn_b_sq;
//	mNoiseCov(X_by, X_by) = pn_b_sq;
//	mNoiseCov(X_bz, X_bz) = pn_b_sq;
//
//	// terrain random walk noise ((m/s)/sqrt(hz)), scales with velocity
//	float pn_t_noise_density =
//		_pn_t_noise_density.get() +
//		(tbl->t_max_grade / 100.0f) * sqrtf(mStateVec(X_vx) * mStateVec(X_vx) + mStateVec(X_vy) * mStateVec(X_vy));
//	mNoiseCov(X_tz, X_tz) = pn_t_noise_density * pn_t_noise_density;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* PE initialization                                              */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int32 PE::InitApp()
{
    int32  iStatus   = CFE_SUCCESS;
    int8   hasEvents = 0;

    iStatus = InitEvent();
    if (iStatus != CFE_SUCCESS)
    {
        (void) CFE_ES_WriteToSysLog("PE - Failed to init events (0x%08lX)\n", iStatus);
        goto PE_InitApp_Exit_Tag;
    }
    else
    {
        hasEvents = 1;
    }

    iStatus = InitPipe();
    if (iStatus != CFE_SUCCESS)
    {
        goto PE_InitApp_Exit_Tag;
    }

    InitData();

    iStatus = InitConfigTbl();
    if (iStatus != CFE_SUCCESS)
    {
        goto PE_InitApp_Exit_Tag;
    }

PE_InitApp_Exit_Tag:
    if (iStatus == CFE_SUCCESS)
    {
        (void) CFE_EVS_SendEvent(PE_INIT_INF_EID, CFE_EVS_INFORMATION,
                                 "Initialized.  Version %d.%d.%d.%d",
								 PE_MAJOR_VERSION,
								 PE_MINOR_VERSION,
								 PE_REVISION,
								 PE_MISSION_REV);
    }
    else
    {
        if (hasEvents == 1)
        {
            (void) CFE_ES_WriteToSysLog("PE - Application failed to initialize\n");
        }
    }

    return iStatus;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* Receive and Process Messages                                    */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int32 PE::RcvSchPipeMsg(int32 iBlocking)
{
    int32           iStatus=CFE_SUCCESS;
    CFE_SB_Msg_t*   MsgPtr=NULL;
    CFE_SB_MsgId_t  MsgId;

    /* Stop Performance Log entry */
    CFE_ES_PerfLogExit(PE_MAIN_TASK_PERF_ID);

    /* Wait for WakeUp messages from scheduler */
    iStatus = CFE_SB_RcvMsg(&MsgPtr, SchPipeId, iBlocking);

    /* Start Performance Log entry */
    CFE_ES_PerfLogEntry(PE_MAIN_TASK_PERF_ID);

    if (iStatus == CFE_SUCCESS)
    {
        MsgId = CFE_SB_GetMsgId(MsgPtr);
        switch (MsgId)
        {
            case PE_WAKEUP_MID:
                /* TODO:  Do something here. */
                break;

            case PE_SEND_HK_MID:
                ProcessCmdPipe();
                ReportHousekeeping();
                break;

            case PX4_VEHICLE_GPS_POSITION_MID:
                memcpy(&mVehicleGpsPositionMsg, MsgPtr, sizeof(mVehicleGpsPositionMsg));
                break;

            case PX4_VEHICLE_STATUS_MID:
                memcpy(&mVehicleStatusMsg, MsgPtr, sizeof(mVehicleStatusMsg));
                break;

            case PX4_VEHICLE_LAND_DETECTED_MID:
                memcpy(&mVehicleLandDetectedMsg, MsgPtr, sizeof(mVehicleLandDetectedMsg));
                break;

            case PX4_ACTUATOR_ARMED_MID:
                memcpy(&mActuatorArmedMsg, MsgPtr, sizeof(mActuatorArmedMsg));
                break;

            case PX4_VEHICLE_ATTITUDE_MID:
                memcpy(&mVehicleAttitudeMsg, MsgPtr, sizeof(mVehicleAttitudeMsg));
                break;

            case PX4_VEHICLE_CONTROL_MODE_MID:
                memcpy(&mVehicleControlModeMsg, MsgPtr, sizeof(mVehicleControlModeMsg));
                break;

            case PX4_SENSOR_COMBINED_MID:
                memcpy(&mSensorCombinedMsg, MsgPtr, sizeof(mSensorCombinedMsg));
                break;

            case PX4_VEHICLE_ATTITUDE_SETPOINT_MID:
                memcpy(&mVehicleAttitudeSetpointMsg, MsgPtr, sizeof(mVehicleAttitudeSetpointMsg));
                break;

            case PX4_MANUAL_CONTROL_SETPOINT_MID:
                memcpy(&mManualControlSetpointMsg, MsgPtr, sizeof(mManualControlSetpointMsg));
                break;

            case PX4_DISTANCE_SENSOR_MID:
                memcpy(&mDistanceSensorMsg, MsgPtr, sizeof(mDistanceSensorMsg));
                break;

            default:
                (void) CFE_EVS_SendEvent(PE_MSGID_ERR_EID, CFE_EVS_ERROR,
                     "Recvd invalid SCH msgId (0x%04X)", MsgId);
        }
    }
    else if (iStatus == CFE_SB_NO_MESSAGE)
    {
        /* TODO: If there's no incoming message, you can do something here, or 
         * nothing.  Note, this section is dead code only if the iBlocking arg
         * is CFE_SB_PEND_FOREVER. */
        iStatus = CFE_SUCCESS;
    }
    else if (iStatus == CFE_SB_TIME_OUT)
    {
        /* TODO: If there's no incoming message within a specified time (via the
         * iBlocking arg, you can do something here, or nothing.  
         * Note, this section is dead code only if the iBlocking arg
         * is CFE_SB_PEND_FOREVER. */
        iStatus = CFE_SUCCESS;
    }
    else
    {
        (void) CFE_EVS_SendEvent(PE_RCVMSG_ERR_EID, CFE_EVS_ERROR,
			  "SCH pipe read error (0x%08lX).", iStatus);
    }

    return iStatus;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* Process Incoming Commands                                       */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void PE::ProcessCmdPipe()
{
    int32 iStatus = CFE_SUCCESS;
    CFE_SB_Msg_t*   CmdMsgPtr=NULL;
    CFE_SB_MsgId_t  CmdMsgId;

    /* Process command messages until the pipe is empty */
    while (1)
    {
        iStatus = CFE_SB_RcvMsg(&CmdMsgPtr, CmdPipeId, CFE_SB_POLL);
        if(iStatus == CFE_SUCCESS)
        {
            CmdMsgId = CFE_SB_GetMsgId(CmdMsgPtr);
            switch (CmdMsgId)
            {
                case PE_CMD_MID:
                    ProcessAppCmds(CmdMsgPtr);
                    break;

                default:
                    /* Bump the command error counter for an unknown command.
                     * (This should only occur if it was subscribed to with this
                     *  pipe, but not handled in this switch-case.) */
                    HkTlm.usCmdErrCnt++;
                    (void) CFE_EVS_SendEvent(PE_MSGID_ERR_EID, CFE_EVS_ERROR,
                                      "Recvd invalid CMD msgId (0x%04X)", (unsigned short)CmdMsgId);
                    break;
            }
        }
        else if (iStatus == CFE_SB_NO_MESSAGE)
        {
            break;
        }
        else
        {
            (void) CFE_EVS_SendEvent(PE_RCVMSG_ERR_EID, CFE_EVS_ERROR,
                  "CMD pipe read error (0x%08lX)", iStatus);
            break;
        }
    }
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* Process PE Commands                                            */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void PE::ProcessAppCmds(CFE_SB_Msg_t* MsgPtr)
{
    uint32  uiCmdCode=0;

    if (MsgPtr != NULL)
    {
        uiCmdCode = CFE_SB_GetCmdCode(MsgPtr);
        switch (uiCmdCode)
        {
            case PE_NOOP_CC:
                HkTlm.usCmdCnt++;
                (void) CFE_EVS_SendEvent(PE_CMD_NOOP_EID, CFE_EVS_INFORMATION,
					"Recvd NOOP. Version %d.%d.%d.%d",
					PE_MAJOR_VERSION,
					PE_MINOR_VERSION,
					PE_REVISION,
					PE_MISSION_REV);
                break;

            case PE_RESET_CC:
                HkTlm.usCmdCnt = 0;
                HkTlm.usCmdErrCnt = 0;
                break;

            default:
                HkTlm.usCmdErrCnt++;
                (void) CFE_EVS_SendEvent(PE_CC_ERR_EID, CFE_EVS_ERROR,
                                  "Recvd invalid command code (%u)", (unsigned int)uiCmdCode);
                break;
        }
    }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* Send PE Housekeeping                                           */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void PE::ReportHousekeeping()
{
    CFE_SB_TimeStampMsg((CFE_SB_Msg_t*)&HkTlm);
    CFE_SB_SendMsg((CFE_SB_Msg_t*)&HkTlm);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* Publish Output Data                                             */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void PE::SendVehicleLocalPositionMsg()
{
    CFE_SB_TimeStampMsg((CFE_SB_Msg_t*)&mVehicleLocalPositionMsg);
    CFE_SB_SendMsg((CFE_SB_Msg_t*)&mVehicleLocalPositionMsg);
}
void PE::SendEstimatorStatusMsg()
{
    CFE_SB_TimeStampMsg((CFE_SB_Msg_t*)&mEstimatorStatusMsg);
    CFE_SB_SendMsg((CFE_SB_Msg_t*)&mEstimatorStatusMsg);
}
void PE::SendVehicleGlobalPositionMsg()
{
    CFE_SB_TimeStampMsg((CFE_SB_Msg_t*)&mVehicleGlobalPositionMsg);
    CFE_SB_SendMsg((CFE_SB_Msg_t*)&mVehicleGlobalPositionMsg);
}
//void PE::SendEkf2InnovationsMsg()
//{
//    CFE_SB_TimeStampMsg((CFE_SB_Msg_t*)&Ekf2InnovationsMsg);
//    CFE_SB_SendMsg((CFE_SB_Msg_t*)&Ekf2InnovationsMsg);
//}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* Verify Command Length                                           */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
boolean PE::VerifyCmdLength(CFE_SB_Msg_t* MsgPtr,
                           uint16 usExpectedLen)
{
    boolean bResult  = TRUE;
    uint16  usMsgLen = 0;

    if (MsgPtr != NULL)
    {
        usMsgLen = CFE_SB_GetTotalMsgLength(MsgPtr);

        if (usExpectedLen != usMsgLen)
        {
            bResult = FALSE;
            CFE_SB_MsgId_t MsgId = CFE_SB_GetMsgId(MsgPtr);
            uint16 usCmdCode = CFE_SB_GetCmdCode(MsgPtr);

            (void) CFE_EVS_SendEvent(PE_MSGLEN_ERR_EID, CFE_EVS_ERROR,
                              "Rcvd invalid msgLen: msgId=0x%08X, cmdCode=%d, "
                              "msgLen=%d, expectedLen=%d",
                              MsgId, usCmdCode, usMsgLen, usExpectedLen);
            HkTlm.usCmdErrCnt++;
        }
    }

    return bResult;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* PE Application C style main entry point.                       */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
extern "C" void PE_AppMain()
{
    oPE.AppMain();
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 */
/* PE Application C++ style main entry point.                     */
/*                                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void PE::AppMain()
{
    /* Register the application with Executive Services */
    uiRunStatus = CFE_ES_APP_RUN;

    int32 iStatus = CFE_ES_RegisterApp();
    if (iStatus != CFE_SUCCESS)
    {
        (void) CFE_ES_WriteToSysLog("PE - Failed to register the app (0x%08lX)\n", iStatus);
    }

    /* Start Performance Log entry */
    CFE_ES_PerfLogEntry(PE_MAIN_TASK_PERF_ID);

    /* Perform application initializations */
    if (iStatus == CFE_SUCCESS)
    {
        iStatus = InitApp();
    }

    if (iStatus == CFE_SUCCESS)
    {
        /* Do not perform performance monitoring on startup sync */
        CFE_ES_PerfLogExit(PE_MAIN_TASK_PERF_ID);
        CFE_ES_WaitForStartupSync(PE_STARTUP_TIMEOUT_MSEC);
        CFE_ES_PerfLogEntry(PE_MAIN_TASK_PERF_ID);
    }
    else
    {
        uiRunStatus = CFE_ES_APP_ERROR;
    }

    /* Application main loop */
    while (CFE_ES_RunLoop(&uiRunStatus) == TRUE)
    {
        RcvSchPipeMsg(PE_SCH_PIPE_PEND_TIME);

        iStatus = AcquireConfigPointers();
        if(iStatus != CFE_SUCCESS)
        {
            /* We apparently tried to load a new table but failed.  Terminate the application. */
            uiRunStatus = CFE_ES_APP_ERROR;
        }
    }

    /* Stop Performance Log entry */
    CFE_ES_PerfLogExit(PE_MAIN_TASK_PERF_ID);

    /* Exit the application */
    CFE_ES_ExitApp(uiRunStatus);
}









/************************/
/*  End of File Comment */
/************************/
