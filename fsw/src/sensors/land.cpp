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
		m_LandTimeout = false;

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
	// measure land
	math::Vector3F y;
	y.Zero();

	if (landMeasure(y) != CFE_SUCCESS)
	{
		return;
	}

	// measurement matrix
	math::Matrix3F10 C;
	C.Zero();
	// y = -(z - tz)
	C[Y_land_vx][X_vx] = 1;
	C[Y_land_vy][X_vy] = 1;
	C[Y_land_agl][X_z] = -1; // measured altitude, negative down dir.
	C[Y_land_agl][X_tz] = 1; // measured altitude, negative down dir.

	OS_printf("C\n");
    for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 10; j++)
		{
			OS_printf("%f, ", C[i][j]);
		}
		OS_printf("\n");
	}
    OS_printf("\n");

	// use parameter covariance
	math::Matrix3F3 R;
	R.Zero();
	R[Y_land_vx][Y_land_vx] = m_Params.LAND_VXY_STDDEV * m_Params.LAND_VXY_STDDEV;
	R[Y_land_vy][Y_land_vy] = m_Params.LAND_VXY_STDDEV * m_Params.LAND_VXY_STDDEV;
	R[Y_land_agl][Y_land_agl] = m_Params.LAND_Z_STDDEV * m_Params.LAND_Z_STDDEV;

	OS_printf("R\n");
    for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			OS_printf("%f, ", R[i][j]);
		}
		OS_printf("\n");
	}
    OS_printf("\n");

	// residual
	//math::Matrix3F3 S_I = inv<float, n_y_land>((C * m_StateCov * C.Transpose()) + R);
	math::Matrix3F3 S_I;
	S_I.Zero();
	S_I = (C * (m_StateCov * C.Transpose())) + R;

	OS_printf("sp\n");
	    for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				OS_printf("%f, ", S_I[i][j]);
			}
			OS_printf("\n");
		}
	OS_printf("\n");

	S_I = S_I.Inversed();

	OS_printf("S_I\n");
    for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			OS_printf("%f, ", S_I[i][j]);
		}
		OS_printf("\n");
	}
    OS_printf("\n");

	math::Vector3F r;
	r.Zero();
	r = C * m_StateVec;

	OS_printf("PRE r\n");
	for(int i = 0; i < 3; i++)
	{
		OS_printf("%f, ", r[i]);
	}
	OS_printf("\n");

	r = y - r;

	OS_printf("r\n");
	for(int i = 0; i < 3; i++)
	{
		OS_printf("%f, ", r[i]);
	}
	OS_printf("\n");

	m_Ekf2InnovationsMsg.HaglInnov = r[Y_land_agl];
	m_Ekf2InnovationsMsg.HaglInnovVar = R[Y_land_agl][Y_land_agl];

	// fault detection
	float beta = (r.Transpose() * (S_I * r)); // TODO: Removed [0][0] index here. This should be a dot product right?

	// artifically increase beta threshhold to prevent fault during landing
	float beta_thresh = 1e2f;

	if (beta / BETA_TABLE[n_y_land] > beta_thresh)
	{
		if (!m_LandFault)
		{
			m_LandFault = true;
			(void) CFE_EVS_SendEvent(PE_LAND_FAULT_ERR_EID, CFE_EVS_ERROR,
									 "Land detector fault, beta %5.2f", double(beta));
		}

		// abort correction
		return;
	}
	else if (m_LandFault)
	{
		m_LandFault = false;
		(void) CFE_EVS_SendEvent(PE_LAND_OK_INF_EID, CFE_EVS_ERROR,
								 "Land detector OK");
	}

	// kalman filter correction always for land detector
	math::Matrix10F3 K;
	K.Zero();
	K = m_StateCov * C.Transpose() * S_I;

	OS_printf("K\n");
    for(int i = 0; i < 10; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			OS_printf("%f, ", K[i][j]);
		}
		OS_printf("\n");
	}
    OS_printf("\n");

	math::Vector10F dx;
	dx.Zero();
	dx = K * r;

	OS_printf("dx\n");
	for(int i = 0; i < 10; i++)
	{
		OS_printf("%f, ", dx[i]);
	}
	OS_printf("\n");

	m_StateVec += dx;
	math::Matrix10F10 res;
	res.Zero();
	res = K * C * m_StateCov;
	OS_printf("res\n");
	res.Print();
	OS_printf("\n");
	m_StateCov -= K * C * m_StateCov;
	OS_printf("STATECOV\n");
	m_StateCov.Print();
	OS_printf("STATE VEC\n");
	for(int i = 0; i < 10; i++)
	{
		OS_printf("%f, ", m_StateVec[i]);
	}
	OS_printf("\n");
}

void PE::landCheckTimeout()
{
	if (m_Timestamp - m_TimeLastLand > LAND_TIMEOUT)
	{
		if (!m_LandTimeout)
		{
			m_LandTimeout = true;
			m_LandCount = 0;
			(void) CFE_EVS_SendEvent(PE_LAND_TIMEOUT_ERR_EID, CFE_EVS_ERROR,
									 "Land detector timeout");
		}
	}
}

