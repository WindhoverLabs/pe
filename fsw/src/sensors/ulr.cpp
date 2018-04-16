#include "../pe_app.h"
#include <math/Matrix1F10.hpp>
#include <math/Matrix1F1.hpp>

/* required number of samples for sensor to initialize */
#define REQ_ULR_INIT_COUNT   (100)
/* 0.1 s */
#define ULR_TIMEOUT          (100000)

void PE::ulrInit()
{
    /* Measure */
	math::Vector1F y;

	if (ulrMeasure(y) != CFE_SUCCESS)
	{
		m_UlrStats.reset();
		return;
	}

	/* If finished */
	if (m_UlrStats.getCount() > REQ_ULR_INIT_COUNT)
	{
		m_UlrAltOrigin = m_UlrStats.getMean()[0];

		(void) CFE_EVS_SendEvent(PE_ULR_OK_INF_EID, CFE_EVS_INFORMATION,
								 "ULR initialized. Mean: (%d) Std dev: (%d) cm",
								 (int)m_UlrStats.getMean()[0],
								 (int)(100 * m_UlrStats.getStdDev()[0]));

		m_UlrTimeout = FALSE;

		if (!m_AltOriginInitialized)
		{
			m_AltOriginInitialized = TRUE;
			m_AltOrigin = m_UlrAltOrigin;
		}
	}
}


int32 PE::ulrMeasure(math::Vector1F &y)
{
	int32 Status = CFE_SUCCESS;
	float d = m_DistanceSensor.CurrentDistance;
	float eps = 0.02f; // 2 cm
	float min_dist = m_DistanceSensor.MinDistance + eps;
	float max_dist = m_DistanceSensor.MaxDistance - eps;

	// prevent driver from setting min dist below eps
	if (min_dist < eps)
	{
		min_dist = eps;
	}
	OS_printf("dist %f\n", d);
	OS_printf("eps %f\n", eps);
	OS_printf("min_dist %f\n", min_dist);
	OS_printf("max_dist %f\n", max_dist);

	// check for bad data
	if (d > max_dist || d < min_dist)
	{
		Status = -1;
		OS_printf("bad dist \n");
		goto ulrMeasure_Exit_Tag;
	}

	// Check for horizontal speed


	/* Measure */
	y.Zero();
	y[0] = (d + m_Params.ULR_OFF_Z) * cosf(m_Euler[0]) * cosf(m_Euler[1]);
	m_UlrStats.update(y);
	m_TimeLastUlr = m_DistanceSensor.Timestamp;

ulrMeasure_Exit_Tag:
	return Status;
}


void PE::ulrCorrect()
{
    CFE_ES_PerfLogEntry(PE_SENSOR_ULR_PERF_ID);
    float cov = 0.0f;

    if (ulrMeasure(m_Ulr.y) != CFE_SUCCESS)
    {
        goto end_of_function;
    }

    /* subtract ulr origin alt */
    m_Ulr.y[0] -= m_UlrAltOrigin;

    /* measured altitude, negative down dir */
    m_Ulr.C[Y_ulr_z][X_z] = -1.0f;
    m_Ulr.C[Y_ulr_z][X_tz] = 1.0f;

    cov = m_DistanceSensor.Covariance;
    if (cov < 1.0e-3f)
    {
    	m_Ulr.R[0][0] = m_Params.ULR_STDDEV * m_Params.ULR_STDDEV;
    }
    else
    {
    	m_Ulr.R[0][0] = cov;
    }

    /* residual */
    /* ((1x10 * 10x10) * 10x1) + 1x1) */
    m_Ulr.S_I = m_Ulr.C * m_StateCov * m_Ulr.C.Transpose() + m_Ulr.R;
    
    /* Take the inverse of a 1x1 matrix (reciprical of the single entry) */
    m_Ulr.S_I[0][0] = 1 / m_Ulr.S_I[0][0];

    /* Vector1F -  (1x10 * Vector10F) */
    m_Ulr.r = m_Ulr.y - (m_Ulr.C * m_StateVec);

    /* fault detection 1F * 1x1 * 1F */
    m_Ulr.beta = m_Ulr.r[0] * m_Ulr.S_I[0][0] * m_Ulr.r[0];

    if (m_Ulr.beta > BETA_TABLE[n_y_ulr])
    {
        if (!m_UlrFault)
        {
            if(Initialized())
            {
                (void) CFE_EVS_SendEvent(PE_ULR_FAULT_ERR_EID, CFE_EVS_ERROR,
                        "Ulr fault, r %5.2f m, beta %5.2f", m_Ulr.r[0], m_Ulr.beta);
            }
            m_UlrFault = TRUE;
        }

    }
    else if (m_UlrFault)
    {
    	m_UlrFault = FALSE;
        (void) CFE_EVS_SendEvent(PE_ULR_OK_INF_EID, CFE_EVS_INFORMATION,
                "Ulr OK");

        m_UlrInitialized = TRUE;
    }

    /* kalman filter correction */

	/* 10x10 * 10x1 * 1x1 */
	m_Ulr.K = m_StateCov * m_Ulr.C.Transpose() * m_Ulr.S_I;

	/* 10x1 * 1x1 */
	m_Ulr.temp = m_Ulr.K * m_Ulr.r;

	m_Ulr.dx = m_Ulr.temp.ToVector();

	/* 10F + 10F*/
	m_StateVec = m_StateVec + m_Ulr.dx;

	/* 10x10 - 10x1 * 1x10 * 10x10 */
	m_StateCov = m_StateCov - m_Ulr.K * m_Ulr.C * m_StateCov;
end_of_function:

    CFE_ES_PerfLogExit(PE_SENSOR_ULR_PERF_ID);
}


void PE::ulrCheckTimeout()
{
    uint64 Timestamp = 0;

	if (m_Timestamp > m_TimeLastUlr)
	{
        Timestamp = m_Timestamp - m_TimeLastUlr;
    }
	else if (m_Timestamp < m_TimeLastUlr)
	{
        Timestamp = m_TimeLastUlr - m_Timestamp;
    }

	if (Timestamp > ULR_TIMEOUT)
	{
		if (!m_UlrTimeout)
		{
			m_UlrTimeout = TRUE;
			m_UlrStats.reset();
			(void) CFE_EVS_SendEvent(PE_ULR_TIMEOUT_ERR_EID, CFE_EVS_ERROR,
									 "Ulr timeout: %llu us", Timestamp);
		}
	}
}
