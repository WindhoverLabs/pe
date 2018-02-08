/************************************************************************
** Includes
*************************************************************************/
#include <string.h>

#include "cfe.h"

#include "pe_app.h"
#include "pe_msg.h"
#include "pe_version.h"

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
	CFE_SB_InitMsg(&m_VehicleLocalPositionMsg,
			PX4_VEHICLE_LOCAL_POSITION_MID, sizeof(PX4_VehicleLocalPositionMsg_t), TRUE);
	CFE_SB_InitMsg(&m_EstimatorStatusMsg,
			PX4_ESTIMATOR_STATUS_MID, sizeof(PX4_EstimatorStatusMsg_t), TRUE);
	CFE_SB_InitMsg(&m_VehicleGlobalPositionMsg,
			PX4_VEHICLE_GLOBAL_POSITION_MID, sizeof(PX4_VehicleGlobalPositionMsg_t), TRUE);
	CFE_SB_InitMsg(&m_Ekf2InnovationsMsg,
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
//	m_BaroStats();
//	m_GpsStats();

	/* Block delays? todo */
	//mXDelay(this, "");
	//mTDelay(this, "");

	/*  */
	memset(&m_Polls, 0, sizeof(m_Polls));
	m_TimeStamp = PX4LIB_GetPX4TimeUs();
	m_TimeStampLastBaro = PX4LIB_GetPX4TimeUs();
	m_TimeLastBaro = 0;
	m_TimeLastGps = 0;
	m_TimeLastLand = 0;

	// reference altitudes
	m_AltOrigin = 0.0;
	m_AltOriginInitialized = FALSE;
	m_BaroAltOrigin = 0.0;
	m_GpsAltOrigin = 0.0;

	// status
	m_ReceivedGps = FALSE;
	m_LastArmedState = FALSE;

	// masks todo
	m_SensorTimeout = 255;
	m_SensorFault = 0;
	m_EstimatorInitialized = 0;

	// state space
	m_StateVec.Zero();
	m_InputVec.Zero();

	InitStateSpace();

	/* Map */
	//mMap_ref.init_done = false;

}

void PE::initStateCov()
{
	m_StateCov.Zero();
	// initialize to twice valid condition
	m_StateCov[X_x][X_x] = 2 * EST_STDDEV_XY_VALID * EST_STDDEV_XY_VALID;
	m_StateCov[X_y][X_y] = 2 * EST_STDDEV_XY_VALID * EST_STDDEV_XY_VALID;
	m_StateCov[X_z][X_z] = 2 * EST_STDDEV_Z_VALID * EST_STDDEV_Z_VALID;
	m_StateCov[X_vx][X_vx] = 2 * ConfigTblPtr->VXY_PUB_THRESH * ConfigTblPtr->VXY_PUB_THRESH;
	m_StateCov[X_vy][X_vy] = 2 * ConfigTblPtr->VXY_PUB_THRESH * ConfigTblPtr->VXY_PUB_THRESH;
	// use vxy thresh for vz init as well
	m_StateCov[X_vz][X_vz] = 2 * ConfigTblPtr->VXY_PUB_THRESH * ConfigTblPtr->VXY_PUB_THRESH;
	// initialize bias uncertainty to small values to keep them stable
	m_StateCov[X_bx][X_bx] = 1e-6;
	m_StateCov[X_by][X_by] = 1e-6;
	m_StateCov[X_bz][X_bz] = 1e-6;
	m_StateCov[X_tz][X_tz] = 2 * EST_STDDEV_TZ_VALID * EST_STDDEV_TZ_VALID;
}

void PE::InitStateSpace()
{
	initStateCov();

	// dynamics matrix
	m_DynamicsMat.Zero();
	// derivative of position is velocity
	m_DynamicsMat[X_x][X_vx] = 1;
	m_DynamicsMat[X_y][X_vy] = 1;
	m_DynamicsMat[X_z][X_vz] = 1;

	// input matrix
	m_InputMat.Zero();
	m_InputMat[X_vx][U_ax] = 1;
	m_InputMat[X_vy][U_ay] = 1;
	m_InputMat[X_vz][U_az] = 1;

	// update components that depend on current state
	updateStateSpace();
	updateStateSpaceParams();
}

void PE::updateStateSpace()
{
	// derivative of velocity is accelerometer acceleration
	// (in input matrix) - bias (in body frame)
	m_DynamicsMat[X_vx][X_bx] = -_R_att[0][0];
	m_DynamicsMat[X_vx][X_by] = -_R_att[0][1];
	m_DynamicsMat[X_vx][X_bz] = -_R_att[0][2];

	m_DynamicsMat[X_vy][X_bx] = -_R_att[1][0];
	m_DynamicsMat[X_vy][X_by] = -_R_att[1][1];
	m_DynamicsMat[X_vy][X_bz] = -_R_att[1][2];

	m_DynamicsMat[X_vz][X_bx] = -_R_att[2][0];
	m_DynamicsMat[X_vz][X_by] = -_R_att[2][1];
	m_DynamicsMat[X_vz][X_bz] = -_R_att[2][2];
}

void PE::updateStateSpaceParams()
{
	// input noise covariance matrix
	m_InputCov.Zero();
	m_InputCov[U_ax][U_ax] = ConfigTblPtr->ACCEL_XY_STDDEV * ConfigTblPtr->ACCEL_XY_STDDEV;
	m_InputCov[U_ay][U_ay] = ConfigTblPtr->ACCEL_XY_STDDEV * ConfigTblPtr->ACCEL_XY_STDDEV;
	m_InputCov[U_az][U_az] = ConfigTblPtr->ACCEL_Z_STDDEV * ConfigTblPtr->ACCEL_Z_STDDEV;

	// process noise power matrix
	m_NoiseCov.Zero();
	float pn_p_sq = ConfigTblPtr->PN_P_NOISE_DENSITY * ConfigTblPtr->PN_P_NOISE_DENSITY;
	float pn_v_sq = ConfigTblPtr->PN_V_NOISE_DENSITY * ConfigTblPtr->PN_V_NOISE_DENSITY;
	m_NoiseCov[X_x][X_x] = pn_p_sq;
	m_NoiseCov[X_y][X_y] = pn_p_sq;
	m_NoiseCov[X_z][X_z] = pn_p_sq;
	m_NoiseCov[X_vx][X_vx] = pn_v_sq;
	m_NoiseCov[X_vy][X_vy] = pn_v_sq;
	m_NoiseCov[X_vz][X_vz] = pn_v_sq;

	// technically, the noise is in the body frame,
	// but the components are all the same, so
	// ignoring for now
	float pn_b_sq = ConfigTblPtr->PN_B_NOISE_DENSITY * ConfigTblPtr->PN_B_NOISE_DENSITY;
	m_NoiseCov[X_bx][X_bx] = pn_b_sq;
	m_NoiseCov[X_by][X_by] = pn_b_sq;
	m_NoiseCov[X_bz][X_bz] = pn_b_sq;

	// terrain random walk noise ((m/s)/sqrt(hz)), scales with velocity
	float pn_t_noise_density =
		ConfigTblPtr->PN_T_NOISE_DENSITY +
		(ConfigTblPtr->T_MAX_GRADE / 100.0f) * sqrtf(m_StateVec[X_vx] * m_StateVec[X_vx] + m_StateVec[X_vy] * m_StateVec[X_vy]);
	m_NoiseCov[X_tz][X_tz] = pn_t_noise_density * pn_t_noise_density;
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
                memcpy(&m_VehicleGpsPositionMsg, MsgPtr, sizeof(m_VehicleGpsPositionMsg));
                break;

            case PX4_VEHICLE_STATUS_MID:
                memcpy(&m_VehicleStatusMsg, MsgPtr, sizeof(m_VehicleStatusMsg));
                break;

            case PX4_VEHICLE_LAND_DETECTED_MID:
                memcpy(&m_VehicleLandDetectedMsg, MsgPtr, sizeof(m_VehicleLandDetectedMsg));
                break;

            case PX4_ACTUATOR_ARMED_MID:
                memcpy(&m_ActuatorArmedMsg, MsgPtr, sizeof(m_ActuatorArmedMsg));
                break;

            case PX4_VEHICLE_ATTITUDE_MID:
                memcpy(&m_VehicleAttitudeMsg, MsgPtr, sizeof(m_VehicleAttitudeMsg));
                break;

            case PX4_VEHICLE_CONTROL_MODE_MID:
                memcpy(&m_VehicleControlModeMsg, MsgPtr, sizeof(m_VehicleControlModeMsg));
                break;

            case PX4_SENSOR_COMBINED_MID:
                memcpy(&m_SensorCombinedMsg, MsgPtr, sizeof(m_SensorCombinedMsg));
                break;

            case PX4_VEHICLE_ATTITUDE_SETPOINT_MID:
                memcpy(&m_VehicleAttitudeSetpointMsg, MsgPtr, sizeof(m_VehicleAttitudeSetpointMsg));
                break;

            case PX4_MANUAL_CONTROL_SETPOINT_MID:
                memcpy(&m_ManualControlSetpointMsg, MsgPtr, sizeof(m_ManualControlSetpointMsg));
                break;

            case PX4_DISTANCE_SENSOR_MID:
                memcpy(&m_DistanceSensorMsg, MsgPtr, sizeof(m_DistanceSensorMsg));
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
    CFE_SB_TimeStampMsg((CFE_SB_Msg_t*)&m_VehicleLocalPositionMsg);
    CFE_SB_SendMsg((CFE_SB_Msg_t*)&m_VehicleLocalPositionMsg);
}
void PE::SendEstimatorStatusMsg()
{
    CFE_SB_TimeStampMsg((CFE_SB_Msg_t*)&m_EstimatorStatusMsg);
    CFE_SB_SendMsg((CFE_SB_Msg_t*)&m_EstimatorStatusMsg);
}
void PE::SendVehicleGlobalPositionMsg()
{
    CFE_SB_TimeStampMsg((CFE_SB_Msg_t*)&m_VehicleGlobalPositionMsg);
    CFE_SB_SendMsg((CFE_SB_Msg_t*)&m_VehicleGlobalPositionMsg);
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
