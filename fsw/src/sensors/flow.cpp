#include "../pe_app.h"
#include <math/Vector2F.hpp>
#include <math/Matrix2F10.hpp>
#include <math/Matrix10F2.hpp>
#include <math/Matrix2F2.hpp>

/* required number of samples for sensor to initialize */
#define REQ_FLOW_INIT_COUNT   (10)
/* 1 s */
#define FLOW_TIMEOUT          (1000000)

#define FLOW_BETA_MAX         (700)
#define FLOW_MIN_AGL         (0.5f) // TODO: Make param
#define FLOW_MIN_QUALITY         (1) // TODO: Make param


void PE::flowInit()
{
    /* Measure */
	math::Vector2F y;

	if (flowMeasure(y) != CFE_SUCCESS)
	{
		m_FlowQStats.reset();
		return;
	}

	/* If finished */
	if (m_FlowStats.getCount() > REQ_FLOW_INIT_COUNT)
	{
		(void) CFE_EVS_SendEvent(PE_FLOW_OK_INF_EID, CFE_EVS_INFORMATION,
								 "Flow initialized. Mean: (%d), Std dev: (%d)",
								 (int)m_FlowQStats.getMean()[0],
								 (int)m_FlowQStats.getStdDev()[0]);

		m_FlowTimeout = FALSE;
	}
}


int32 PE::flowMeasure(math::Vector2F &y)
{
	int32 Status = CFE_SUCCESS;
	math::Euler euler;
	float d = 0;
	float flow_x_rad = 0;
	float flow_y_rad = 0;
	float dt_flow = 0;
	float gyro_x_rad = 0;
	float gyro_y_rad = 0;
	math::Vector3F delta_b;
	math::Vector3F delta_n;

	/* Check for sane pitch/roll */
	if (m_Euler[0] > 0.5f || m_Euler[1] > 0.5f) {
		Status = -1;
		goto flowMeasure_Exit_Tag;
	}

	// check for agl
	if (m_AglLowPass.m_State < FLOW_MIN_AGL) {
		Status = -1;
		goto flowMeasure_Exit_Tag;
	}

	// check quality
	if (m_OpticalFlowMsg.Quality < FLOW_MIN_QUALITY) {
		Status = -1;
		goto flowMeasure_Exit_Tag;
	}

	// calculate range to center of image for flow
	if (!m_TzEstValid) {
		Status = -1;
		goto flowMeasure_Exit_Tag;
	}

	//euler = (math::Quaternion(m_VehicleAttitudeMsg.Q)).ToEuler(); TODO

	d = m_AglLowPass.m_State * cosf(euler[0]) * cosf(euler[1]);

	// optical flow in x, y axis
	// TODO consider making flow scale a state of the kalman filter
	flow_x_rad = m_OpticalFlowMsg.PixelFlowXIntegral;// TODO * _flow_scale.get();
	flow_y_rad = m_OpticalFlowMsg.PixelFlowYIntegral;//TODO * _flow_scale.get();
	dt_flow = m_OpticalFlowMsg.IntegrationTimespan / 1.0e6f;

	if (dt_flow > 0.5f || dt_flow < 1.0e-6f) {
		Status = -1;
		goto flowMeasure_Exit_Tag;
	}

	// angular rotation in x, y axis
	gyro_x_rad = 0;
	gyro_y_rad = 0;

//	if (_fusion.get() & FUSE_FLOW_GYRO_COMP) {
//		gyro_x_rad = _flow_gyro_x_high_pass.update(
//					 _sub_flow.get().gyro_x_rate_integral);
//		gyro_y_rad = _flow_gyro_y_high_pass.update(
//					 _sub_flow.get().gyro_y_rate_integral);
//	}

	// compute velocities in body frame using ground distance
	// note that the integral rates in the optical_flow uORB topic are RH rotations about body axes
	delta_b[0] = float(fabs(flow_y_rad - gyro_y_rad) * d);
	delta_b[1] = -(flow_x_rad - gyro_x_rad) * d;
	delta_b[2] = 0.0f;

	// rotation of flow from body to nav frame
	delta_n = m_RotationMat * delta_b;

	// imporant to timestamp flow even if distance is bad
	m_TimeLastFlow = m_Timestamp;

	// measurement
	y[Y_flow_vx] = delta_n[0] / dt_flow;
	y[Y_flow_vy] = delta_n[1] / dt_flow;

	//m_FlowQStats.update(Scalarf(m_OpticalFlowMsg.Quality));
	m_FlowQStats.update(float(m_OpticalFlowMsg.Quality));

flowMeasure_Exit_Tag:
	return Status;
}


void PE::flowCorrect()
{
    CFE_ES_PerfLogEntry(PE_SENSOR_FLOW_PERF_ID);
    const float h_min = 2.0f;
	const float h_max = 8.0f;
	const float v_min = 0.5f;
	const float v_max = 1.0f;
	const float p[5] = {0.04005232f, -0.00656446f, -0.26265873f,  0.13686658f, -0.00397357f};
    float h = 0;
    float v = 0;
    float flow_vxy_stddev = 0;
    float rotrate_sq = 0;
    float rot_sq = 0;

    if (flowMeasure(m_Flow.y) != CFE_SUCCESS)
    {
        goto end_of_function;
    }

    /* measured altitude, negative down dir */
    m_Flow.C[Y_flow_vx][X_vx] = 1.0f;
    m_Flow.C[Y_flow_vy][X_vy] = 1.0f;

    // polynomial noise model, found using least squares fit
	// h, h**2, v, v*h, v*h**2

	// prevent extrapolation past end of polynomial fit by bounding independent variables
	h = m_AglLowPass.m_State;
	v = sqrtf(m_Flow.y * m_Flow.y); //was norm()

	if (h > h_max) {
		h = h_max;
	}

	if (h < h_min) {
		h = h_min;
	}

	if (v > v_max) {
		v = v_max;
	}

	if (v < v_min) {
		v = v_min;
	}

	// compute polynomial value
	flow_vxy_stddev = p[0] * h + p[1] * h * h + p[2] * v + p[3] * v * h + p[4] * v * h * h;

	rotrate_sq = m_VehicleAttitudeMsg.RollSpeed * m_VehicleAttitudeMsg.RollSpeed
			   + m_VehicleAttitudeMsg.PitchSpeed * m_VehicleAttitudeMsg.PitchSpeed
			   + m_VehicleAttitudeMsg.YawSpeed * m_VehicleAttitudeMsg.YawSpeed;

	rot_sq = m_Euler[0] * m_Euler[0] + m_Euler[1] * m_Euler[1];

    /* Vector2F -  (1x10 * Vector10F) */
    m_Flow.r = m_Flow.y - (m_Flow.C * m_StateVec);

    /* residual */
    /* ((2x10 * 10x10) * 10x2) + 2x2) */
    m_Flow.S_I = m_Flow.C * m_StateCov * m_Flow.C.Transpose() + m_Flow.R;
    m_Flow.S_I = m_Flow.S_I.Inversed();

    /* fault detection 1F * 1x1 * 1F */
    m_Flow.beta = (m_Flow.r.Transpose() * (m_Flow.S_I * m_Flow.r))[0][0];

    if (m_Flow.beta > FLOW_BETA_MAX)
    {
        if (!m_FlowFault)
        {
            if(Initialized())
            {
                (void) CFE_EVS_SendEvent(PE_FLOW_FAULT_ERR_EID, CFE_EVS_ERROR,
                        "Flow fault, beta %5.2f", m_Flow.beta);
            }
            m_FlowFault = TRUE;
        }
        goto end_of_function;
    }
    else if (m_FlowFault)
    {
    	m_FlowFault = FALSE;
        (void) CFE_EVS_SendEvent(PE_FLOW_OK_INF_EID, CFE_EVS_INFORMATION,
                "Flow OK, beta %5.2f", m_Flow.beta);

        m_FlowInitialized = TRUE;
    }

    /* kalman filter correction */
    if (!m_FlowFault)
    {
		/* 10x10 * 10x2 * 2x2 */
		m_Flow.K = m_StateCov * m_Flow.C.Transpose() * m_Flow.S_I;

		/* 10x2 * 2x2 */
		m_Flow.temp = m_Flow.K * m_Flow.r;

		m_Flow.dx = m_Flow.temp.ToVector();

		/* 10F + 10F*/
		m_StateVec = m_StateVec + m_Flow.dx;

		/* 10x10 - 10x2 * 2x10 * 10x10 */
		m_StateCov = m_StateCov - m_Flow.K * m_Flow.C * m_StateCov;
    }

end_of_function:
    CFE_ES_PerfLogExit(PE_SENSOR_FLOW_PERF_ID);
}


void PE::flowCheckTimeout()
{
    uint64 Timestamp = 0;

	if (m_Timestamp > m_TimeLastFlow)
	{
        Timestamp = m_Timestamp - m_TimeLastFlow;
    }
	else if (m_Timestamp < m_TimeLastFlow)
	{
        Timestamp = m_TimeLastFlow - m_Timestamp;
    }

	if (Timestamp > FLOW_TIMEOUT)
	{
		if (!m_FlowTimeout)
		{
			m_FlowTimeout = TRUE;
			m_FlowStats.reset();
			(void) CFE_EVS_SendEvent(PE_FLOW_TIMEOUT_ERR_EID, CFE_EVS_ERROR,
									 "Flow timeout: %llu us", Timestamp);
		}
	}
}
