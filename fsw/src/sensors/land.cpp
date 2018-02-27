#include "../pe_app.h"

// required number of samples for sensor
// to initialize
//
#define 	REQ_LAND_INIT_COUNT (1)
#define 	LAND_TIMEOUT   		(1000000) // 1.0 s

void PE::landInit()
{
	// measure
	math::Vector3F y;
	y.Zero();

	if (landMeasure(y) != CFE_SUCCESS)
	{
		m_LandCount = 0;
	}

	// if finished
	if (m_LandCount > REQ_LAND_INIT_COUNT)
	{
		m_LandTimeout = FALSE;

		(void) CFE_EVS_SendEvent(PE_LAND_OK_INF_EID, CFE_EVS_INFORMATION,
								 "Land detector initialized");
	}
}

int PE::landMeasure(math::Vector3F &y)
{
	m_TimeLastLand = m_Timestamp;
	y.Zero();
	m_LandCount += 1;
	return CFE_SUCCESS;
}

void PE::landCorrect()
{
    CFE_ES_PerfLogEntry(PE_SENSOR_LAND_PERF_ID);
	// measure land
	//math::Vector3F y;
	//y.Zero();

	if (landMeasure(m_Land.y) != CFE_SUCCESS)
	{
		goto end_of_function;
	}

	// measurement matrix
	//math::Matrix3F10 C;
	//C.Zero();
	// y = -(z - tz)
	m_Land.C[Y_land_vx][X_vx] = 1.0f;
	m_Land.C[Y_land_vy][X_vy] = 1.0f;
	m_Land.C[Y_land_agl][X_z] = -1.0f; // measured altitude, negative down dir.
	m_Land.C[Y_land_agl][X_tz] = 1.0f; // measured altitude, negative down dir.

	// use parameter covariance
	//math::Matrix3F3 R;
	//R.Zero();
	m_Land.R[Y_land_vx][Y_land_vx] = m_Params.LAND_VXY_STDDEV * m_Params.LAND_VXY_STDDEV;
	m_Land.R[Y_land_vy][Y_land_vy] = m_Params.LAND_VXY_STDDEV * m_Params.LAND_VXY_STDDEV;
	m_Land.R[Y_land_agl][Y_land_agl] = m_Params.LAND_Z_STDDEV * m_Params.LAND_Z_STDDEV;

	// residual
	//math::Matrix3F3 S_I;
	//S_I.Zero();
	m_Land.S_I = (m_Land.C * (m_StateCov * m_Land.C.Transpose())) + m_Land.R;
	m_Land.S_I = m_Land.S_I.Inversed();

	//math::Vector3F r;
	//r.Zero();
	m_Land.r = m_Land.y - m_Land.C * m_StateVec;

	m_Ekf2InnovationsMsg.HaglInnov = m_Land.r[Y_land_agl];
	m_Ekf2InnovationsMsg.HaglInnovVar = m_Land.R[Y_land_agl][Y_land_agl];

	// fault detection
	m_Land.beta = (m_Land.r.Transpose() * (m_Land.S_I * m_Land.r));

	// artifically increase beta threshhold to prevent fault during landing
	m_Land.beta_thresh = 1e2f;

	if (m_Land.beta / BETA_TABLE[n_y_land] > m_Land.beta_thresh)
	{
		if (!m_LandFault)
		{
			m_LandFault = TRUE;
			(void) CFE_EVS_SendEvent(PE_LAND_FAULT_ERR_EID, CFE_EVS_ERROR,
									 "Land detector fault, beta %5.2f", double(m_Land.beta));
		}

		// abort correction
		goto end_of_function;
	}
	else if (m_LandFault)
	{
		m_LandFault = FALSE;
		(void) CFE_EVS_SendEvent(PE_LAND_OK_INF_EID, CFE_EVS_ERROR,
								 "Land detector OK");
		m_LandInitialized = TRUE;
	}

	// kalman filter correction always for land detector
	//math::Matrix10F3 K;
	//K.Zero();
	m_Land.K = m_StateCov * m_Land.C.Transpose() * m_Land.S_I;

	//math::Vector10F dx;
	//dx.Zero();
	m_Land.dx = m_Land.K * m_Land.r;

	m_StateVec = m_StateVec + m_Land.dx;
    //OS_printf("PRE LAND\n");
    //m_StateCov.Print();
	m_StateCov = m_StateCov - m_Land.K * m_Land.C * m_StateCov;
    //OS_printf("LAND CORRECTED\n");
    //m_StateCov.Print();
end_of_function:

    CFE_ES_PerfLogExit(PE_SENSOR_LAND_PERF_ID);
}

void PE::landCheckTimeout()
{
	if (m_Timestamp - m_TimeLastLand > LAND_TIMEOUT)
	{
		if (!m_LandTimeout)
		{
			m_LandTimeout = TRUE;
			m_LandCount = 0;
			(void) CFE_EVS_SendEvent(PE_LAND_TIMEOUT_ERR_EID, CFE_EVS_ERROR,
									 "Land detector timeout");
		}
	}
}

